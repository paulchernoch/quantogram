use std::cell::RefCell;
use std::ops::Bound::{Included,Excluded};
use std::fmt::{Formatter,Debug,Result};
use float_extras::f64::frexp;
use std::convert::Into;
use skiplist::SkipMap;
mod half_sample_mode;
use half_sample_mode::{HalfSampleModeCache};


/// A bin to be used inside a Histogram with a midpoint value and a weight.
#[derive(Copy, Clone, Debug)]
struct HistogramBin {
    /// Midpoint value for the bin.
    pub sample : f64,

    /// Sample weight, which must be positive.
    /// If zero, it indicates a removed bin.
    pub weight : f64
}

impl HistogramBin {
    /// Construct a new HistogramBin.
    pub fn new(sample: f64, weight: f64) -> Self {
        HistogramBin {
            sample: sample,
            weight: weight
        }
    }

    /// Validates that the weight is non-negative and finite, else function will panic.
    fn validate(&self) {
        if self.weight < 0.0 || !self.weight.is_finite() {
            panic!("HistogramBin weight must be a non-negative, finite number.");
        }
    }

    /// Copy values from another HistogramBin into this one.
    pub fn set(&mut self, copy_from: Self) {
        self.sample = copy_from.sample;
        self.weight = copy_from.weight;
        self.validate();
    }
}

/// Categorize numbers as negative, zero, positive, NAN or an infinity
enum SampleCategory {
    Negative,
    Zero,
    Positive,
    NegativeInfinity,
    PositiveInfinity,
    NotANumber
}

impl SampleCategory {
    /// Decide which category of number a sample is.
    pub fn categorize(sample : f64) -> Self {
        match sample {
            // These arms are sorted by most to least common which improved performance noticeably. 
            x if x > 0.0 && x != f64::INFINITY => SampleCategory::Positive,
            x if x < 0.0 && x != f64::NEG_INFINITY => SampleCategory::Negative,
            x if x == 0.0 => SampleCategory::Zero,
            x if x == f64::INFINITY => SampleCategory::PositiveInfinity,
            x if x == f64::NEG_INFINITY => SampleCategory::NegativeInfinity,
            _ => SampleCategory::NotANumber
        }
    }
}

/// Provides a weighted Histogram of f64 values for computing approximate quantiles.
/// This guarantees a configurable maximum absolute relative error
/// and uses sparse storage to reduce memory usage.
/// 
/// Worst case accuracy defaults to one percent (0.01) absolute relative error.
/// The error is unbiased, uniform for the entire range of numbers.
/// The error for quantiles 0 and 1 (the minimum and maximum, respectively) is guaranteed to 
/// be zero, except if either of those values is removed.
/// 
/// If all inserted values are given a weight of one,
/// this behaves as an unweighted (normal) histogram. 
/// 
/// Samples may be added or removed. 
/// However, removing a sample that equals the minimum or maximum will cause
/// those values to be replaced by the center value of the appropriate extreme histogram bucket.
/// 
/// The full valid range of floats is divided into two levels of buckets.
/// For the default case with 1% error, here is what that means:
/// 
///    Top. The top level divides the number range into buckets for each power of two and between
///         positive and negative numbers. 
///         It has a maximum of 508 buckets in two sparse lists. 
///         f64 values can range from 2^-126 to 2^127, thus there are a maximum of 
///         254 buckets for positive numbers and 254 buckets for negatives. 
/// Bottom. The second level of buckets are in a single sparse list spanning the full range 
///         from the smallest negative number to the largest positive number. 
///         Each power of two range is broken into 35 buckets whose size varies exponentially. 
///         Each bucket is larger than the previous by a factor of 1.02
///         (or smaller by 1.02, over the range of negative values). 
///         The value of 1.02 was chosen because 1.02^35 = 1.999889553, which differs 
///         from 2 by 0.00011. That means that the 35th bucket for each power of two is slightly larger than
///         the rest, but not by enough to wreck the error guarantee.
///         There are a maximum of 1 + 508*35 buckets in this list, or 17,781. 
///         (The "1" is for the zero bucket.) 
/// 
/// A bucket is not created until at least one value is added to it.
/// Removing the last item in a bucket will not cause its memory to be freed; 
/// its weight will be set to zero. 
/// 
/// (If you are familiar with the Rule of 70 used to estimate the doubling period 
/// for a given interest rate, it was probably chosen because of this nice property of 1.02.)
/// 
/// The error rate of 0.01 and the bin scale factor of 1.02 are related in this way: 
/// ```text
/// 
///             (1 + error)     (1 + 1/101)      102
///    scale = ------------- = ------------- = ------- = 1.02
///             (1 - error)     (1 - 1/101)      100
/// ```
/// So technically, the worst case error is 1/101, or 0.99%, not 1%, but the math is easier
/// when using 1.02 instead of 1.020202020202 and the error on the last bucket is ideal. 
/// (A smaller error in the last bucket is available for error = 1/176, but the memory requirements
/// are much higher.)
/// 
/// Usage: 
/// 
///   Typical usage (unweighted samples with the default accuracy of 1%): 
/// 
///      use quantogram::Quantogram;
///      let mut q = Quantogram::new();
///      q.add(10.0);
///      q.add(40.0);
///      q.add(20.0);
///      q.add(30.0);
///      q.add(50.0);
/// 
///      assert_eq!(q.min().unwrap(), 10.0);
///      assert_eq!(q.max().unwrap(), 50.0);
///      assert_eq!(q.mean().unwrap(), 30.0);
///      assert_eq!(q.median().unwrap(), 30.0);
///      assert_eq!(q.quantile(0.75).unwrap(), 40.0);
/// 
///      q.remove(10.0);
///      q.remove(20.0);
///      assert_eq!(q.mean().unwrap(), 40.0);
///   
/// 
/// Notes: 
///   1. Coarse bins are for powers of two not some other value because 
///      getting the largest power of two less than or equal to a number calls a fast intrinsic function.
///      This makes assigning a number to a bin very fast.
/// 
///   2. When inquiring about a quantile, a value will be returned so long as one
///      number in range is added. NANs and Infinities will be disregarded. 
///      The NANs and Infinities are available as separate counts.  
/// 
///   3. Unbounded errors are possible in the edge case of a large gap in the middle of the data.
///      Take the case of the median. If there are an even number of items in the data with a 
///      large gap in the exact middle, then the proper formula for the median is the mean 
///      of the last value below the gap and the first value above the gap.
///      To correct for this, use the fussy_quantile method. 
///      It will probe for quantiles at Φ + ε and Φ - ε for a small value of ε, 
///      and average the pair that best span the gap. 
///      If the histogram is unweighted (all weights are one) then the value of
///      ε should be 1/2N, where N is the number of items already added to the Histogram.
///      If the samples are weighted, not sure what to do.
/// 
///   4. If one sample is in a bin, the value for the bin will be set to the accurate sample,
///      not the midpoint of the range for that bin. If a second sample is added, the
///      bin value will be set to the midpoint. Consequently, bins with one
///      sample added have no error. 
///      
pub struct Quantogram {
    /// Exponential Scaling factor that relates the start value of one fine bin to the next.
    /// value[N+1] = value[N] * growth
    /// The default is 1.02, which yields an absolute relative error of 1%.
    /// This number must be in the range (1,√2]. 
    /// Smaller values guarantee smaller errors but require more storage space and execution time.
    growth : f64,

    /// Number of bins used to store values per power of two range of numbers.
    /// This must equal log(2, growth), rounded down.
    /// For example, if growth is 1.02, this must be 35.
    /// The larger this value, the more memory is required.
    bins_per_doubling : usize,

    /// Memoize powi(growth, n) for n in [0,bins_per_doubling)
    growth_bin_power_memo : Vec<f64>,

    /// Unweighted count of all insertions minus all removals.
    total_count : usize,
     
    /// Total weight of all inserted values other than NAN or +/- Infinity. 
    total_weight : f64,

    /// Weighted sum of all samples (excluding NANs and infinities).
    /// This permits the exact weighted mean to be calculated.
    weighted_sum : f64,

    /// Running weighted variance
    running_variance : f64,

    /// Total weight of all inserted NAN values. 
    nan_weight : f64,

    /// Total weight of all inserted +Infinity values. 
    plus_ininity_weight : f64,

    /// Total weight of all inserted -Infinity values. 
    minus_ininity_weight : f64,

    /// Total weight of all negative values (included in total_weight).
    negative_weight : f64,

    /// Total weight of all zeroes (included in total_weight).
    zero_weight : f64,
    
    /// Total weight of all positive values (included in total_weight).
    positive_weight : f64,
   
    /// Maximum added sample that is not +Infinity, initialized to NAN. 
    maximum : f64,

    /// Minimum added sample that is not -Infinity, initialized to NAN. 
    minimum : f64,

    /// This remains true until a finite number that is not a number is added.
    /// When true, this causes quantile to round its result.
    only_integers : bool,

    /// To reduce memory usage, force all numbers whose magnitude is
    /// greater than this power of two to be counted as +Infinity (overflow). 
    /// Default is -(2^127), the highest supported by IEEE-754 64 bit floats. 
    negative_overflow : f64,

    /// To reduce memory usage, force all numbers whose magnitude
    /// is less than this power of two to be counted as zero (underflow). 
    /// Default is -(2^-126), the lowest supported by IEEE-754 64 bit floats. 
    negative_underflow : f64,

    /// To reduce memory usage, force all numbers whose magnitude is
    /// greater than this power of two to be counted as +Infinity (overflow). 
    /// Default is +(2^127), the highest supported by IEEE-754 64 bit floats. 
    positive_overflow : f64,
 
    /// To reduce memory usage, force all numbers whose magnitude
    /// is less than this power of two to be counted as zero (underflow). 
    /// Default is +(2^-126), the lowest supported by IEEE-754 64 bit floats. 
    positive_underflow : f64,

    /// Bins for powers of two of positive numbers.
    /// The key is the power of two, from -126 to 127.
    /// 
    /// Example: If the key is 2, the bin is for all numbers in the range [4,8).
    positive_coarse_bins : SkipMap<isize, HistogramBin>,

    /// Bins for powers of two of negative numbers.
    /// The key is the power of two, from -126 to 127.
    /// 
    /// Example: If the key is 3, the bin is for all numbers in the range (-16,-8].
    negative_coarse_bins : SkipMap<isize, HistogramBin>,

    /// Bins for zeroes, positive and negative numbers.
    /// The key is related to the power of two and the position within the 35 bins per
    /// power of two range: 
    ///     key = ((power + 126)*35 + bin + 1)*sign
    /// where:
    ///     power = the power of two, from -126 to 127
    ///     bin   = 0 to 34
    ///     sign  = 1 for positive numbers, -1 for negative numbers. 
    /// Zeroes are always stored at key = 0.
    fine_bins : SkipMap<isize, HistogramBin>,

    /// Natural Log of the grown factor, precomputed for speed.
    log_of_growth : f64,

    /// Cache of mode candidates for computing the mode.
    mode_cache : ModeCache,

    /// Cache for use with hsm (half-sample mode). 
    hsm_cache : RefCell<HalfSampleModeCache>
}

impl Quantogram {

    // ////////////////////////////
    //                           //
    //       Constructors        //
    //                           //
    // ////////////////////////////

    /// Create a Quantogram with the default error rate of 1%.
    pub fn new() -> Self {
        let bins = 35;
        let growth = 1.02;
        let mut q = Quantogram {
            growth : growth,
            bins_per_doubling : bins,

            growth_bin_power_memo : Vec::new(),

            total_count : 0,
            total_weight : 0.0,
            weighted_sum : 0.0,
            running_variance : 0.0,
            nan_weight : 0.0, 
            plus_ininity_weight : 0.0,
            minus_ininity_weight : 0.0,
            negative_weight : 0.0,
            zero_weight : 0.0,
            positive_weight : 0.0,
            maximum : f64::NAN,
            minimum : f64::NAN, 

            only_integers : true,

            negative_underflow : -(2.0_f64.powi(-126)),
            negative_overflow : -(2.0_f64.powi(127)),
            positive_underflow : 2.0_f64.powi(-126),
            positive_overflow : 2.0_f64.powi(127),

            // The SkipMaps do not preallocate bins. This lets it know
            // how many levels to make each skiplist to make it optimal.
            positive_coarse_bins : SkipMap::with_capacity(254),
            negative_coarse_bins : SkipMap::with_capacity(254),
            fine_bins : SkipMap::with_capacity(17641),
            log_of_growth : growth.ln(),
            mode_cache : ModeCache::new(),
            hsm_cache : RefCell::new(HalfSampleModeCache::new_default())
        };
        q.memoize_growth_bin_power();
        q
    } 

    /// Create a Quantogram where growth and bins are set, and underflow and overflow bounds are set
    /// to the given powers of two. 
    ///   - Inserted values having absolute magnitude below 2^smallest_power will be treated as zeroes.
    ///   - Inserted values having absolute magnitude above 2^largest_power will be treated as infinities.
    /// 
    /// Note: To ensure that consist values for growth and bins are set, it is advised to use QuantogramBuilder 
    ///       instead of directly calling this constructor.
    pub fn with_configuration(growth: f64, bins: usize, smallest_power: isize, largest_power: isize) -> Self {
        let valid_power = -126..=127;
        if !valid_power.contains(&smallest_power) || !valid_power.contains(&largest_power) {
            panic!("Power of two constraints must be in range [-126,127]");
        }
        if smallest_power >= largest_power {
            panic!("Power of two constraints not given in ascending order");
        }
        let mut q = Quantogram {
            growth : growth,
            bins_per_doubling : bins,

            growth_bin_power_memo : Vec::new(),

            total_count : 0,
            total_weight : 0.0,
            weighted_sum : 0.0,
            running_variance : 0.0,
            nan_weight : 0.0, 
            plus_ininity_weight : 0.0,
            minus_ininity_weight : 0.0,
            negative_weight : 0.0,
            zero_weight : 0.0,
            positive_weight : 0.0,
            maximum : f64::NAN,
            minimum : f64::NAN, 

            only_integers : true,

            negative_underflow : -(2.0_f64.powi(smallest_power as i32)),
            negative_overflow : -(2.0_f64.powi(largest_power as i32)),
            positive_underflow : 2.0_f64.powi(smallest_power as i32),
            positive_overflow : 2.0_f64.powi(largest_power as i32),

            // The SkipMaps do not preallocate bins. This lets it know
            // how many levels to make each skiplist to make it optimal.
            positive_coarse_bins : SkipMap::with_capacity(254),
            negative_coarse_bins : SkipMap::with_capacity(254),
            // If the range is constrained, fewer bins will be required,
            // so the capacity can be reduced. 
            fine_bins : SkipMap::with_capacity(((largest_power - smallest_power + 1)*35 + 1) as usize),
            log_of_growth : growth.ln(),
            mode_cache : ModeCache::new(),
            hsm_cache : RefCell::new(HalfSampleModeCache::new_default())
        };
        q.memoize_growth_bin_power();
        q
    }

    pub fn replace_hsm_cache(&mut self, new_cache: HalfSampleModeCache) {
        self.hsm_cache = RefCell::new(new_cache);
        self.hsm_cache.borrow_mut().clear();
    }

    // ////////////////////////////
    //                           //
    //    Add/Remove Samples     //
    //                           //
    // ////////////////////////////

    /// Add a sample to the histogram with a weight of one.
    pub fn add(&mut self, sample: f64) {
        self.add_weighted(sample, 1.0);
    }

    /// Remove a sample from the histogram, deducting one from the weight.
    /// If the weight goes negative, the call panics.
    pub fn remove(&mut self, sample: f64) {
        self.add_weighted(sample, -1.0);
    }
    
    /// Add a sample to the histogram with the given weight.
    /// If the weight is positive, the sample is added, otherwise removed.
    /// If the cumulative weight for the bin holding that sample goes negative, 
    /// the call panics.
    pub fn add_weighted(&mut self, sample: f64, weight: f64) -> f64 {
        self.hsm_cache.borrow_mut().record(sample);
        let bounded_sample = self.over_or_underflow(sample);
        let adjusted_fine_weight;
        self.adjust_aggregates(bounded_sample, weight);
        match SampleCategory::categorize(bounded_sample) {
            SampleCategory::Positive => {
                self.positive_weight += weight;
                let (fine_key, fine_low, fine_midpoint, fine_high) = self.get_fine_key_with_midpoint(bounded_sample, self.growth, self.bins_per_doubling).unwrap();
                adjusted_fine_weight = Self::increment_key_weight(&mut self.fine_bins, fine_key, fine_midpoint, sample, weight); 
                self.adjust_mode(fine_key, fine_low, fine_midpoint, fine_high, adjusted_fine_weight);
                let (coarse_key, coarse_midpoint) = Self::get_coarse_key_with_midpoint(bounded_sample).unwrap();
                Self::increment_key_weight(&mut self.positive_coarse_bins, coarse_key, coarse_midpoint, sample, weight);
            },
            SampleCategory::Zero => {
                self.zero_weight += weight;
                let (fine_key, fine_low, fine_midpoint, fine_high) = self.get_fine_key_with_midpoint(bounded_sample, self.growth, self.bins_per_doubling).unwrap();
                adjusted_fine_weight = Self::increment_key_weight(&mut self.fine_bins, fine_key, fine_midpoint, sample, weight); 
                self.adjust_mode(fine_key, fine_low, fine_midpoint, fine_high, adjusted_fine_weight);
            },
            SampleCategory::Negative => {
                self.negative_weight += weight;
                let (fine_key, fine_low, fine_midpoint, fine_high) = self.get_fine_key_with_midpoint(bounded_sample, self.growth, self.bins_per_doubling).unwrap();
                adjusted_fine_weight = Self::increment_key_weight(&mut self.fine_bins, fine_key, fine_midpoint, sample, weight); 
                self.adjust_mode(fine_key, fine_low, fine_midpoint, fine_high, adjusted_fine_weight);
                let (coarse_key, coarse_midpoint) = Self::get_coarse_key_with_midpoint(bounded_sample).unwrap();
                Self::increment_key_weight(&mut self.negative_coarse_bins, coarse_key, coarse_midpoint, sample, weight);
            },
            SampleCategory::NotANumber => {
                self.nan_weight += weight;
                adjusted_fine_weight = self.nan_weight;
            },
            SampleCategory::PositiveInfinity => {
                self.plus_ininity_weight += weight;
                adjusted_fine_weight = self.plus_ininity_weight;
            },
            SampleCategory::NegativeInfinity => {
                self.minus_ininity_weight += weight;
                adjusted_fine_weight = self.minus_ininity_weight;
            }
        }
        adjusted_fine_weight
    }

    /// Add many samples to the Quantogram, all having a weight of 1.0.
    pub fn add_unweighted_samples<'a, S>(&mut self, samples: impl Iterator<Item = &'a S>) 
    where S: 'a + Into<f64> + Copy
    {
        for sample in samples {
            self.add((*sample).into());
        }
    }

    // ////////////////////////////
    //                           //
    //      Basic Measures       //
    //                           //
    //   mean min max count      //
    //   median variance stddev  //
    //   mode hsm                //
    //   quantile quantile_at    //
    //   fussy_quantile          //
    //   finite zero nan         //
    //                           //
    // ////////////////////////////


    /// Return actual (not estimated) weighted mean of inserted samples (with machine precision),
    /// or None if no finite values have been inserted.
    /// Excludes all NANs and Infinities that have been inserted.
    pub fn mean(&self) -> Option<f64> {
        if self.total_weight > 0.0 { Some(self.weighted_sum / self.total_weight) }
        else { None }
    }

    /// Return actual (not estimated) minimum of inserted samples (with machine precision),
    /// or None if no finite values have been inserted.
    /// Excludes all NANs and Infinities that have been inserted.
    pub fn min(&self) -> Option<f64> {
        if self.minimum.is_finite() { Some(self.minimum) }
        else { None }
    }

    /// Return actual (not estimated) maximum of inserted samples (with machine precision),
    /// or None if no finite values have been inserted.
    /// Excludes all NANs and Infinities that have been inserted.
    pub fn max(&self) -> Option<f64> {
        if self.maximum.is_finite() { Some(self.maximum) }
        else { None }
    }

    /// Return count of inserted samples, including NANs and Infinities.
    pub fn count(&self) -> usize {
        self.total_count
    }

    /// Return the weighted fraction of values that are finite (not NAN or +/- Infinity)
    /// as a value in the range [0,1]. 
    /// If no values have yet been added, return NAN. 
    pub fn finite(&self) -> f64 {
        if self.total_count == 0 {
            f64::NAN
        }
        else {
            let numerator = self.total_weight;
            let denominator = self.added_weight();
            numerator / denominator
        }
    }

    /// Return the weighted fraction of values that are zero, including NANs and infinities in the basis. 
    /// If no values have yet been added, return NAN. 
    pub fn zero(&self) -> f64 {
        if self.total_count == 0 {
            f64::NAN
        }
        else {
            let numerator = self.zero_weight;
            let denominator = self.added_weight();
            numerator / denominator
        }
    }

    /// Return the weighted fraction of values that are NAN, including infinities in the basis. 
    /// If no values have yet been added, return NAN. 
    pub fn nan(&self) -> f64 {
        if self.total_count == 0 {
            f64::NAN
        }
        else {
            let numerator = self.nan_weight;
            let denominator = self.added_weight();
            numerator / denominator
        }
    }

    /// Return an estimate of the median.
    pub fn median(&self) -> Option<f64> {
        self.quantile(0.5)
    }

    /// Exact weighted variance of all added finite values.
    pub fn variance(&self) -> f64 {
        self.running_variance
    }

    /// Exact weighted population standard deviation of all added finite values.
    pub fn stddev(&self) -> Option<f64> {
        if self.total_weight <= 0.0 {
            None
        }
        else {
            Some((self.running_variance / self.total_weight).sqrt())
        }
    }

    /// Return an estimate of the mode.
    /// If no finite samples have been added, returns an empty list. 
    /// If there are multiple modes tieing with the same weight,
    /// return all of them in a list.
    /// If all samples in the collection are integers, round the mode.
    pub fn mode(&self) -> Vec<f64> {
        let unrounded_mode = self.mode_cache.mode(&self.fine_bins);
        if self.only_integers {
            let rounded_mode = unrounded_mode
                .into_iter()
                .map(|x| x.round())
                .collect();
            rounded_mode
        }
        else {
            unrounded_mode
        }
    }

    /// Estimate the half-sample mode.
    /// More resistant than the mode to outliers and contamination (noise). 
    /// 
    /// Based on the mode estimator by Robertson and Cryer (1974), and an 
    /// algorithm described in "On a Fast, Robust Estimator of the Mode" by David Bickel and Rudolf Fruhwirth.
    /// That algorithm is applied to raw samples, whereas this is applied 
    /// to already histogrammed values.
    /// 
    /// Edge case: If fewer than five bins of data are present in the histogram,
    /// the true mode will be returned, unless the data is multi-moded, in which case
    /// None will be returned.
    /// 
    /// NOTE: 
    pub fn hsm(&self) -> Option<f64> {
        if self.total_weight == 0.0 { None }
        else {
            let mut weights: Vec<f64> = Vec::new();
            let mut samples: Vec<f64> = Vec::new();
            for bin in self.fine_bins.values() {
                samples.push(bin.sample);
                weights.push(bin.weight);
            }

            // Cache may be updated because of interior mutability.
            self.hsm_cache.borrow_mut().hsm(&samples, &weights, self.total_weight)
            
            // If we wanted to bypass the cache, do it this way:
            // half_sample_mode(&samples, &weights, self.total_weight)
        }
    }

    /// Estimate the quantile for the inserted data.
    ///   For the minimum, use phi = 0.0.
    ///   For the median, use phi = 0.5. 
    ///   For the 90th percentile, phi = 0.9.
    ///   For the maximum, use phi = 1.0.
    /// 
    /// Added samples that are NANs or Infinities are excluded from the computation.
    pub fn quantile(&self, phi: f64) -> Option<f64> {
        if phi < 0.0 || phi > 1.0 || self.total_weight <= 0.0 {
            None
        }
        else if phi == 0.0 {
            self.min()
        }
        else if phi == 1.0 {
            self.max()
        }
        else {
            let mut target_cume_weight = self.total_weight * phi;

            // Narrow down the position of the quantile to negative, zero or positive numbers. 
            if target_cume_weight <= self.negative_weight {
                // Quantile falls in the negative range
                // Search two maps: coarse and fine. 
                let (coarse_index, remainder) = Self::find_containing_bucket(target_cume_weight, &self.negative_coarse_bins).unwrap();
                // Tricky, because in the negative range we must go from the highest magnitude to the lowest magnitude.
                // Coarse index is therefore the negative of the power of two only for negative numbers.
                let high_fine_index = ((-coarse_index + 126) * (self.bins_per_doubling as isize) + 0 + 1) * -1 + 1; // Exclusive upper bound
                let low_fine_index = high_fine_index - (self.bins_per_doubling as isize) + 1; // Inclusive lower bound
                let (_fine_index, _remainder2, sample_at_quantile) = Self::find_containing_bucket_in_range(low_fine_index, high_fine_index, remainder, &self.fine_bins).unwrap();
                // println!("Negatives. coarse = {}. fine range = [{},{}) Fine found = {}", coarse_index, low_fine_index, high_fine_index, fine_index);
                // Previously required a second call to get, which was slower.
                // let sample_at_quantile = self.fine_bins.get(&fine_index).unwrap().sample;
                return Some(self.conditional_round(sample_at_quantile));
            }
            target_cume_weight -= self.negative_weight;
            if target_cume_weight <= self.zero_weight {
                return Some(0.0);
            }
            else {
                target_cume_weight -= self.zero_weight;
            }

            // Quantile falls in the positive range
            let (coarse_index, remainder) = Self::find_containing_bucket(target_cume_weight, &self.positive_coarse_bins).unwrap();
            let low_fine_index = ((coarse_index + 126) * (self.bins_per_doubling as isize) + 0 + 1) * 1;
            let high_fine_index = low_fine_index + (self.bins_per_doubling as isize); // Exclusive upper bound
            let (_fine_index, _remainder2, sample_at_quantile) 
                = Self::find_containing_bucket_in_range(low_fine_index, high_fine_index, remainder, &self.fine_bins)
                  .unwrap();
            return Some(self.conditional_round(sample_at_quantile));
        }
    }

    /// Quantile measurement that corrects for issues that occur when there is a large gap in samples 
    /// right where the quantile falls.
    /// 
    /// For example, the median of { 1,2,3,4,5,95,96,97,98,99 } is 50! This is because there are an even number of items
    /// in the set, so you have to average the middle two values.
    /// 
    /// Three quantile calculations will be taken, at Φ - ε, Φ and Φ + ε.
    /// If the difference between the results at the two lower quantiles
    /// differs substantially from the difference between the results at the two higher quantiles,
    /// assume there is a gap and average the middle value and the most extreme value. 
    /// 
    /// threshold_ratio decides whether a gap exists. If none exists, no averaging occurs. 
    /// It is used to compare the deltas between quantiles computed at the desired phi and values of 
    /// phi slightly lower or higher. If the deltas are similar, no gap exists. 
    /// A value over 2.0 seems sensible, but testing should be done to determine the best value.
    pub fn fussy_quantile(&self, phi: f64, threshold_ratio: f64) -> Option<f64> {
        if phi < 0.0 || phi > 1.0 || self.total_weight <= 0.0 {
            None
        }
        else if phi == 0.0 {
            self.min()
        }
        else if phi == 1.0 {
            self.max()
        }
        else {
            // TODO: This derivation of epsilon is appropriate for unweighted samples (all having weight = 1.0). 
            // Unsure what the proper value is for weighted samples.
            let epsilon = 1.0 / (2.0 * self.total_count as f64);
            
            let q_middle = self.quantile(phi);
            if phi <= 1.5 * epsilon || phi >= 1.0 - 1.5 * epsilon || q_middle.is_none() {
                q_middle
            } 
            else {
                let q_low = self.quantile(phi - epsilon).unwrap();
                let q_high = self.quantile(phi + epsilon).unwrap();
                let q_median = q_middle.unwrap();

                let low_difference = q_median - q_low;
                let high_difference = q_high - q_median;
                if low_difference >= high_difference * threshold_ratio {
                    Some((q_low + q_median)/2.0)
                }
                else if high_difference >= low_difference * threshold_ratio {
                    Some((q_high + q_median)/2.0)
                }
                else {
                    None
                }     
            }
        }
    }

    /// Get the range of quantiles between which the given value falls,
    /// or None if no finite samples have been added yet or the 
    /// queried value is not finite.
    pub fn quantile_at(&self, value: f64) -> Option<(f64,f64)> {
        let min_opt = self.min();
        if min_opt.is_none() || !value.is_finite() {
            return None;
        } 
        let min = min_opt.unwrap();
        let max = self.max().unwrap();
        match self.count() {
            0 => None,
            1 => {
                if value < min { Some((0.0, 0.0)) }
                else if value > min { Some((1.0, 1.0)) }
                else { Some((0.0, 1.0)) }
            },
            _ if value < min => Some((0.0, 0.0)),
            _ if value > max => Some((1.0, 1.0)),
            _ if min == max => Some((0.0, 1.0)),
            _ if value == min => Some((0.0, 0.0)),
            _ if value == max => Some((1.0, 1.0)),
            _ if value == 0.0 => {
                let neg = self.negative_weight;
                let sum = self.total_weight;
                let zero = self.zero_weight;
                Some((neg / sum, (neg + zero)/sum))
            },
            _ => {
                // TODO: It is possible to use the coarse skipmaps to narrow in on the location
                // in the fine skipmap faster, hence speedup this calculation. 
                // This is a less commonly used method, so optimization is not a priority. 
                let (end_key, _fine_low, _fine_midpoint, _fine_high) 
                    = self.get_fine_key_with_midpoint(value, self.growth, self.bins_per_doubling).unwrap();
                let (start_key, start_weight) = 
                    if value >= 0.0 { (1_isize, self.negative_weight + self.zero_weight) } 
                    else { (isize::MIN , 0.0) };
                let mut cume_weight = start_weight;
                let sum = self.total_weight;
                for (key, bucket) in self.fine_bins.range(Included(&start_key), Included(&isize::MAX)) {
                    let weight = bucket.weight;
                    if *key == end_key {
                        // Value falls within an existing bin, so we must define a range
                        // or phi values.
                        return Some((cume_weight / sum, (cume_weight + weight)/sum));
                    }
                    else if *key > end_key {
                        // Value would go in a bin not yet present in the collection, so we passed by it.
                        // Since value falls between bins, it gets a well defined phi value,
                        // where the bottom and top of range are identical.
                        return Some((cume_weight / sum, cume_weight/sum));
                    }
                    cume_weight += weight;
                }
                // Should never fall through, because the maximum sample always has a bin defined for it, 
                // and we already tested if value > max.
                None
            }
        }
    }

    // ////////////////////////////
    //                           //
    //   Measures of Dispersion  //
    //                           //
    //   range q1 q3 iqr         //
    //   quartile_deviation      //
    //   coeff_of_range          //
    //   coeff_of_quartile_dev   //
    //   coeff_of_stddev         //
    //   coeff_of_variation      //
    //                           //
    // ////////////////////////////

    /// Range between the highest and lowest values sampled.
    pub fn range(&self) -> Option<f64> {
        match (self.max(), self.min()) {
            (Some(max), Some(min)) => Some(max - min),
            _ => None
        }
    }

    /// Sample at the first quartile.
    pub fn q1(&self) -> Option<f64> { self.quantile(0.25) }

    /// Sample at the third quartile.
    pub fn q3(&self) -> Option<f64> { self.quantile(0.75) }

    /// Inter quartile range, which is (Q3 - Q1)
    pub fn iqr(&self) -> Option<f64> {
        match (self.q3(), self.q1()) {
            (Some(q3), Some(q1)) => Some(q3 - q1),
            _ => None
        }  
    }

    /// Quartile deviation or semi-interquartile range, which is (Q3 - Q1) / 2
    pub fn quartile_deviation(&self) -> Option<f64> {
        match (self.q3(), self.q1()) {
            (Some(q3), Some(q1)) => Some((q3 - q1) / 2.0),
            _ => None
        }  
    }

    /// Coefficient of range. 
    /// 
    /// ```text
    ///                        max - min
    ///    coeff of range =  ------------
    ///                        max + min
    /// ```
    pub fn coeff_of_range(&self) -> Option<f64> {
        match (self.max(), self.min()) {
            (Some(max), Some(min)) => {
                let sum = max + min;
                if sum == 0.0 { None }
                else { Some((max - min) / sum) }
            },
            _ => None
        }
    }

    /// Coefficient of quartile deviation, which measures the relative spread between 
    /// the first and third quartiles.
    /// 
    /// ```text
    ///                                    q3 - q1
    ///    coeff of quartile deviation =  ---------
    ///                                    q3 + q1
    /// ```    
    pub fn coeff_of_quartile_dev(&self) -> Option<f64> {
        match (self.q3(), self.q1()) {
            (Some(q3), Some(q1)) => {
                let sum = q3 + q1;
                if sum == 0.0 { None }
                else { Some((q3 - q1) / sum) }
            },
            _ => None
        }
    }

    /// Coefficient of standard deviation.
    /// 
    /// This give the standard deviation divided by the mean.
    pub fn coeff_of_stddev(&self) -> Option<f64> {
        match (self.stddev(), self.mean()) {
            (Some(stddev), Some(mean)) => {
                if mean == 0.0 { None }
                else { Some(stddev / mean) }
            },
            _ => None
        }
    }

    /// Coefficient of variation.
    /// 
    /// This give the standard deviation divided by the mean expressed as a percentage.
    pub fn coeff_of_variation(&self) -> Option<f64> {
        match self.coeff_of_stddev() {
            Some(coeff) => Some(coeff * 100.0),
            _ => None
        }
    }

    // ////////////////////////////
    //                           //
    //       Diagnostics         //
    //                           //
    // ////////////////////////////

    /// Computes a relative measure of the storage being used by the histogram. 
    /// 
    /// As bins are dynamically added to the sparse Skipmaps, this increases. 
    pub fn size(&self) -> usize {
        7 // counts of NANs, +/- Infinity, zeroes, three empty Skipmaps are always kept
        + self.positive_coarse_bins.len()
        + self.negative_coarse_bins.len()
        + self.fine_bins.len()
    }

    // ////////////////////////////
    //                           //
    //    Utility functions      //
    //                           //
    // ////////////////////////////

    /// Compute growth^N for N = 0..bins_per_doubling and store in a Vec
    /// to memoize the result as an optimization, to avoid having to call powi function. 
    fn memoize_growth_bin_power(&mut self) {
        for bin in 0..=self.bins_per_doubling {
            self.growth_bin_power_memo.push(self.growth.powi(bin as i32));
        }
    }

    /// Approximate the log(x, growth) using a Padé Approximant
    /// and take the floor, yielding an integer. 
    /// Good accuracy if x is in the range [1.0,2.0]. (Accuracy actually good up to e = 2.71828.)
    /// 
    /// ```text
    /// A good second-order Padé approximant to ln(x) for the region x ∈ [1,e] is:
    /// 
    ///                3x² - 3
    ///    ln(x) =  -------------
    ///              x² + 4x + 1
    /// 
    /// A third-order Padé approximant to ln(x) for the region x ∈ [1,e] is:
    /// 
    ///              11x³ + 27x² - 27x - 11
    ///    ln(x) =  ------------------------
    ///               3x³ + 27x² + 27x + 3
    /// ```
    /// To transform to use the growth factor as base, we divide by ln(growth). 
    /// 
    /// See https://math.stackexchange.com/questions/2238880/good-approximation-to-lnx-for-x-in-1-x-e
    /// for how the 2nd order Pade approximant was derived.
    /// 
    /// The third-order approximant came from Wolfram Alpha with some algebraic rearrangement 
    /// to optimize the number of arithmetic operations: 
    ///   > PadeApproximant(ln[x], {x, 1, {3, 3}})
    /// 
    /// Over the range [1,2], the worst case absolute relativer error for the:
    ///   2nd order Padé approximant is 0.12%. 
    ///   3nd order Padé approximant is 0.004%. (0.033% for e = 2.71828)
    /// The 3rd order accuracy was needed. In release mode, it is twice as fast as the system log. 
    fn pade_approximate_log_floor(&self, x: f64) -> isize {
        let x2 = x * x;
        let x3 = x2 * x;

        // 2nd order: 
        // let numerator = 3.0 * (x2 - 1.0);
        // let denominator = self.log_of_growth * (x2 + 4.0*x + 1.0);

        // 3rd order (after dividing all terms by 27 and changing base): 
        let numerator = (11.0/27.0) * x3 + x2 - x - (11.0/27.0);
        let denominator = self.log_of_growth * (x3/9.0 + x2 + x + (1.0/9.0));

        let ln_x = numerator / denominator;
        ln_x.floor() as isize
    }

    /// Find the key to the containing bucket by Iterating over buckets in a map 
    /// in increasing order and accumulating the weights until one is reached where the total
    /// falls within the range belonging to the bucket.
    /// Also return the amount of weight remaining prior to the final bucket.
    fn find_containing_bucket(cume_weight: f64, map: &SkipMap<isize, HistogramBin>) -> Option<(isize,f64)> {
        let mut remaining_weight = cume_weight;
        for (k, bucket) in map.iter() {
            let weight = bucket.weight;
            if weight >= remaining_weight {
                return Some((*k, remaining_weight));
            }
            remaining_weight -= weight;
        }
        None
    }

    /// Find the key to the containing bucket in the given key range
    /// by Iterating over buckets in a map 
    /// accumulating the weights until one is reached where the total
    /// falls within the range belonging to the bucket.
    /// 
    /// Returns a tuple with the key, remaining weight, and the bucket sample value.
    fn find_containing_bucket_in_range(start_key: isize, end_key: isize, cume_weight: f64, map: &SkipMap<isize, HistogramBin>) -> Option<(isize,f64,f64)> {
        let mut remaining_weight = cume_weight;
        for (k, bucket) in map.range(Included(&start_key), Excluded(&end_key)) {
            let weight = bucket.weight;
            if weight >= remaining_weight {
                return Some((*k, remaining_weight, bucket.sample));
            }
            remaining_weight -= weight;
        }
        None
    }

    /// Bookkeeping done when adding a new sample or removing an existing one. 
    fn adjust_aggregates(&mut self, sample: f64, weight: f64) {
        self.adjust_minimum(sample);
        self.adjust_maximum(sample);
        let weight_is_finite = weight.is_finite();
        let sample_is_finite = sample.is_finite();
        if weight < 0.0 && weight_is_finite {
            // Negative weights are interpreted as a removal.
            self.total_count -= 1;
        }
        else {
            self.total_count += 1;
        }
        if sample_is_finite && weight_is_finite {
            let previous_mean = self.mean().unwrap_or_default();
            self.total_weight += weight;
            self.weighted_sum += weight * sample;
            let current_mean = self.mean().unwrap_or_default();
            self.running_variance += weight * (sample - previous_mean) * (sample - current_mean);
        }

        if sample_is_finite && sample.fract() != 0.0 {
            // Once we add numbers with fractions, assume that the quantile can have a fraction.
            self.only_integers = false;
        }
    }

    /// Restrict a value based on preset overflow and underflow values. 
    /// Values that underflow are set to zero.
    /// Values that overflow are set to +/- infinity.
    fn over_or_underflow(&self, sample: f64) -> f64 {
        if sample > self.negative_underflow && sample < self.positive_underflow {
            0.0
        }
        else if sample > self.positive_overflow {
            f64::INFINITY
        }
        else if sample < self.negative_overflow {
            f64::NEG_INFINITY
        }
        else {
            sample
        }
    }

    /// If the new value is finite and less than the minimum 
    /// (or the minimum is unset, hence a NAN), update the minimum.
    fn adjust_minimum(&mut self, sample: f64) {
        if self.minimum.is_nan() {
            if sample.is_finite() {
                self.minimum = sample;
            }
        }
        else if sample < self.minimum && sample.is_finite() {
            self.minimum = sample;
        }
    }

    /// If the new value is finite and greater than the maximum
    /// (or the maximum is unset, hence a NAN), update the maximum.
    fn adjust_maximum(&mut self, sample: f64) {
        if self.maximum.is_nan() {
            if sample.is_finite() {
                self.maximum = sample;
            }
        }
        else if sample > self.maximum && sample.is_finite() {
            self.maximum = sample;
        }
    }

    fn adjust_mode(
        &mut self, 
        bin_key: isize, 
        bin_low_value: f64, 
        bin_midpoint: f64, 
        bin_high_value: f64, 
        bin_weight: f64) {
        // TODO: Add a configuration parameter here and in QuantogramBuilder
        // to enable/disable mode calculations. 
        // It is a drag on performance and not all callers need the mode.

        // NOTE: Do not update the mode if the weight is 1. 
        // Otherwise if all values are distinct (common in unit tests)
        // then you store a huge list of modes. This kills performance, too.
        if bin_midpoint.is_finite() && bin_weight > 0.0 && bin_weight != 1.0 {
            self.mode_cache.update(bin_key, bin_low_value, bin_high_value, bin_weight);
        }
    }

    /// Update bin weight for given key, or create new bin for key if none exists.
    /// Returns the updated weight.
    fn increment_key_weight(map: &mut SkipMap<isize, HistogramBin>, key: isize, bucket_sample: f64, accurate_sample: f64, added_weight: f64) -> f64 {
        match map.get_mut(&key) {
            Some(bucket) => {
                // When updating an existing bucket, use the bucket_sample value, which is the midpoint value,
                // unless the previous bucket used the same value as the accurate value.
                let new_bucket:HistogramBin;
                if bucket.sample == accurate_sample {
                    new_bucket = HistogramBin::new(accurate_sample, bucket.weight + added_weight);
                    bucket.set(new_bucket); 
                }
                else {
                    new_bucket = HistogramBin::new(bucket_sample, bucket.weight + added_weight);
                    bucket.set(new_bucket);                    
                }
                new_bucket.weight
            },
            None => {
                // When creating a new bucket, use the accurate_sample. 
                // This means that bins with a single sample are 100% accurate,
                // but when the second sample is added, shift to the midpoint for the bin.
                map.insert(key, HistogramBin::new(accurate_sample, added_weight));
                added_weight
            }
        }
    }

    /// If only integers have been added to the histogram, round the value.
    fn conditional_round(&self, sample: f64) -> f64 {
        if self.only_integers && sample.is_finite() {
            sample.round()
        }
        else {
            sample
        }
    }

    /// Sum of weights of all samples added, including NAN and infinities.
    fn added_weight(&self) -> f64 {
        self.total_weight + self.nan_weight + self.plus_ininity_weight + self.minus_ininity_weight
    }

    
    /// Exponent on the largest power of two less than or equal to the magnitude of the given number,
    /// or None if a NAN, Infinity or Zero.
    /// 
    /// ```
    ///   use quantogram::Quantogram;
    ///   assert_eq!(Quantogram::power_of_two(5.0).unwrap(), 2);
    ///   assert_eq!(Quantogram::power_of_two(4.0).unwrap(), 2);
    ///   assert_eq!(Quantogram::power_of_two(1.5).unwrap(), 0);
    ///   assert_eq!(Quantogram::power_of_two(-8.5).unwrap(), 3);
    ///   assert_eq!(Quantogram::power_of_two(0.0), None);
    ///   assert_eq!(Quantogram::power_of_two(0.25).unwrap(), -2);
    /// ```
    pub fn power_of_two(sample: f64) -> Option<isize> {
        if sample.is_finite() && sample != 0.0 {
            let (_mantissa, exponent) = frexp(sample);
            Some(exponent - 1)
        }
        else {
            None
        }
    }

    /// Break a sample into three parts: 
    ///   - sign is -1 for negative numbers or +1 for positive numbers
    ///   - power is exponent on the largest power of two less than or equal to the absolute value of the number.
    ///   - remainder is the mantissa multiplied by two, which is a number in the semi-open range [1.0,2.0)
    /// 
    /// Returns None if the sample is not finite or is zero.
    fn sign_power_remainder(sample: f64) -> Option<(isize, isize, f64)> {
        if sample.is_finite() && sample != 0.0 {
            let (mantissa, exponent) = frexp(sample);
            if sample > 0.0 {
                Some((1, exponent - 1, mantissa * 2.0))
            }
            else {
                Some((-1, exponent - 1, mantissa.abs() * 2.0))
            }
        }
        else {
            None
        }
    }

    /// Raise two to the given power.
    #[inline(always)]
    fn two_power(power: isize) -> f64 {
        match power {
            0..=63 => (1_usize << power) as f64,
            64..=127 => (1_u128 << power) as f64,
            -63..=-1 => 1.0 / ((1_usize << -power) as f64),
            -126..=-64 => 1.0 / ((1_u128 << -power) as f64),
            _ => f64::NAN
        }
    }

    /// Get the key that maps a sample to its coarse bin. 
    /// For positive numbers, this is the floor of the log base 2 of the sample.
    /// For negative numbers, negate the power of two, so that negative numbers
    /// sort from low to high. (A negative number with a higher power of two is "smaller"
    /// or more negative than one with a lower power of two.)
    /// 
    /// Returns a Tuple with the bin key and the midpoint value for the bin.
    fn get_coarse_key_with_midpoint(sample: f64) -> Option<(isize,f64)> {
        match Self::sign_power_remainder(sample) {
            Some((sign, power, _remainder)) => {
                let low = Self::two_power(power); // 2.0_f64.powi(power as i32);
                let high = low * 2.0;
                let bucket_middle = (sign as f64) * (low + high) / 2.0;
                // Multiply by the sign because negative numbers with large exponents sort before those with small exponents
                Some((sign * power, bucket_middle))
            },
            _ => None
        }
    }

    /// Get the key that maps a number to its fine histogram bin given the desired exponential growth factor.
    /// If the growth factor is 1.02, there will be 35 bins per power of two.
    /// Thus bins_per_doubling must match the growth_factor: 
    /// 
    /// ```text
    ///     growth^bins_per_doubling ~ 2.0.
    /// ```
    /// 
    /// Bins are defined for zeroes, and finite positive and negative numbers.
    /// The key is related to the power of two and the position within the 35 bins per
    /// power of two range: 
    /// ```text
    ///     key = ((power + 126)*35 + bin + 1)*sign
    /// ```
    /// where:
    /// ```text
    ///     power = the power of two, from -126 to 127
    ///     bin   = 0 to 34 (for growth = 1.02)
    ///     sign  = 1 for positive numbers, -1 for negative numbers. 
    /// Zeroes are always stored at key = 0.
    /// ```
    /// Returns a Tuple with the bin key and the low, mid, and high point values for the bin. 
    fn get_fine_key_with_midpoint(&self, sample: f64, growth: f64, bins_per_doubling: usize) -> Option<(isize,f64,f64,f64)> {
        match Self::sign_power_remainder(sample) {
            Some((sign, power, remainder)) => {
                // Use of system log is more accurate but slower. 
                // This inaccuracy adds a slight bias towards zero. 
                // Rare samples will be put into the next lower bin. Since that bin is slightly narrower, 
                // this may actually reduce the error slightly as the estimate is pushed towards a closer bin center. 
                // let bin = remainder.log(growth).floor() as isize; 
                let bin = self.pade_approximate_log_floor(remainder);

                let key = ((power + 126) * (bins_per_doubling as isize) + bin + 1) * sign;
                let two_power = Self::two_power(power); // f64::powi(2.0, power as i32);
                let bucket_low = two_power * self.growth_bin_power_memo[bin as usize]; // Memo of: f64::powi(growth, bin as i32);
                let bucket_high = bucket_low * growth;
                let bucket_middle = (sign as f64) * (bucket_low + bucket_high) / 2.0;
                Some((key, bucket_low, bucket_middle, bucket_high))
            },
            _ => { 
                if sample == 0.0 {
                    Some((0, 0.0, 0.0, 0.0))
                }
                else {
                    None
                }
            }
        }
    }

}

impl Debug for Quantogram {
    fn fmt(&self, f: &mut Formatter) -> Result {
        let mut fine = "(".to_owned();
        for (k, bucket) in self.fine_bins.iter() {
            fine = fine + &format!("\n              {}->{:.4}:{}", k, bucket.sample, bucket.weight);   
        }
        fine = fine + "\n            )";

        let mut coarse_plus = "(".to_owned();
        for (k, bucket) in self.positive_coarse_bins.iter() {
            coarse_plus = coarse_plus + &format!("\n              {}->{:.4}:{}", k, bucket.sample, bucket.weight);   
        }
        coarse_plus = coarse_plus + "\n            )";

        let mut coarse_minus = "(".to_owned();
        for (k, bucket) in self.negative_coarse_bins.iter() {
            coarse_minus = coarse_minus + &format!("\n              {}->{:.4}:{}", k, bucket.sample, bucket.weight);   
        }
        coarse_minus = coarse_minus + "\n            )";

        write!(f, "Quantogram.
          growth = {}, bins per doubling = {}
          total count = {}, weight = {}, sum = {}, variance = {:.3},
          NAN = {}, -Inf = {}, +Inf = {}, Integers? = {}
          - = {}, 0 = {}, + = {}
          min = {}, mean = {:?}, max = {}
          under/overflow = {}/{}
          coarse bins = 
            Minus {}
            Plus {}
          fine bins = {}
        ", 
            self.growth, self.bins_per_doubling,
            self.total_count, self.total_weight, self.weighted_sum,
            self.running_variance,
            self.nan_weight, self.minus_ininity_weight, self.plus_ininity_weight, self.only_integers,
            self.negative_weight, self.zero_weight, self.positive_weight,
            self.minimum, self.mean(), self.maximum,
            self.positive_underflow.log2(), self.positive_overflow.log2(),
            coarse_minus,
            coarse_plus,
            fine
        )
    }
}

/// A fluent API Builder for Quantograms. 
/// 
/// If error, growth or bins_per_doubling are changed, the other two 
/// are modified to be consistent with it. If error or growth imply a non-integral
/// value for bins_per_doubling, take the ceiling of the derived bins_per_ceiling
/// and rederive the other two from that, to satisfy a minimum guarantee. 
/// 
/// No values may cause bins_per_doubling to exceed 1000, hence the error 
/// may not be set lower than 0.0347%.
/// 
/// Example of creating a Quantogram with a maximum of 2% absolute relative error: 
/// 
/// 
///   use quantogram::QuantogramBuilder;
///   let q = QuantogramBuilder().with_error(0.02).build();
///
/// 
/// This is the number of bins used for various error rates: 
/// ```text
///   Bins    Error Rate
///   ----    ----------
///     2        17.2 %   
///     3        11.5 %   
///     4         8.6 %   
///     5         6.9 %   
///     7         4.9 %   
///    12         2.9 %   
///    18         1.9 %   
///    35         1.0 %   
///    70         0.5 %
///   139         0.25%
///   347         0.10%
///   694         0.05%
///   867         0.04%
///  1000         0.0347%
/// ```
#[derive(Copy, Clone, Debug)]
pub struct QuantogramBuilder {
    /// Desired absolute relative error rate.
    /// The valid range is a fraction in the range [0.00035,0.17].
    /// Defaults to approximately 0.01.
    /// Altering error changes growth and bins_per_doubling to be consistent.
    /// 
    /// This is related to the number of bins per doubling (b) 
    /// or growth via these formulas:
    /// ```text             
    ///              1/b
    ///            2      - 1     growth - 1
    ///   error = ------------ = ------------
    ///              1/b          growth + 1
    ///            2      + 1
    /// ```
    error : f64,

    /// Number of bins used to store values per power of two range of numbers.
    /// The valid range is [2,1000].
    /// The larger this value, the more memory is required.
    /// Defaults to 35 and must exceed one.
    /// Altering this changes error and growth to be consistent.
    /// 
    /// The number of bins (b) is related to the growth factor 
    /// via this formula: 
    /// ```text
    ///                1/b
    ///     growth = 2
    /// ```
    bins_per_doubling : usize,

    /// Exponential scale factor relating start value of one bin to the next. 
    /// The valid range is [1.00069,1.4143].
    /// Defaults to approximately 1.02. 
    /// Altering growth changes error and bins_per_doubling to be consistent.
    /// ```text
    ///           b    
    ///     growth  = 2
    /// ```
    growth : f64,

    /// A sample whose magnitude is less than this power of two will underflow and be considered a zero.
    /// Must be in the range [-126,127] but also be less than largest_power.
    /// Defaults to -126.
    smallest_power: isize, 
    
    /// A sample whose magnitude is more than this power of two will overflow and be considered +/- Infinity.
    /// Must be in the range [-126,127] but also be more than smallest_power.
    /// Defaults to +127.
    largest_power: isize,

    hsm_cache: Option<HalfSampleModeCache>
}

impl QuantogramBuilder {
    /// Create a new builder that defaults to an error rate of 1% 
    /// with 35 bins per doubling and a growth factor of 1.02.
    pub fn new() -> Self {
        QuantogramBuilder {
            error: 0.01,
            bins_per_doubling: 35,
            growth: 1.02,
            smallest_power: -126,
            largest_power: 127,
            hsm_cache: None
        }
    }

    /// Build a Quantogram using the collected configuration values.
    pub fn build(self) -> Quantogram {
        let mut q = Quantogram::with_configuration(self.growth, self.bins_per_doubling, self.smallest_power, self.largest_power);
        match self.hsm_cache {
            Some(cache) => { q.replace_hsm_cache(cache); },
            None => (),
        }
        q
    }

    pub fn with_hsm_cache(mut self, hsm_cache: &HalfSampleModeCache)-> Self {
        let mut cleared_cache = hsm_cache.clone();
        cleared_cache.clear();
        self.hsm_cache = Some(cleared_cache);
        self
    } 

    /// Configure the underflow of samples.
    /// A sample whose magnitude has a power of two less than the given values will be set to zero.
    /// 
    /// Example: If you only need to compute quantiles to two decimal places:
    /// ```text
    ///          set power = -7, since 2^-7 = 1/128 < 0.01.
    /// ```
    pub fn with_smallest_power(mut self, power: isize)-> Self {
        self.smallest_power = 
            if power < -126 { -126 }
            else if power >= 127 { 127 }
            else { power };
        self
    }

    /// Configure the overflow of samples.
    /// A sample whose magnitude has a power of two greater than the given values will be set to +/- Infinity.
    /// 
    /// Example: If you only need to study numbers below a thousand, 
    ///          set power = 10, since 2^10 = 1024 > 1000.
    pub fn with_largest_power(mut self, power: isize)-> Self {
        self.largest_power = 
            if power <= -126 { -125 }
            else if power > 127 { 127 }
            else { power };
        self
    }

    /// Configure to target the given maximum absolute relative error.  
    pub fn with_error(mut self, error: f64) -> Self {
        (&mut self).set_from_error(error);
        self
    }

    /// Configure to target the given growth factor from one bucket to the next. 
    /// This will be adjusted to the nearest value that is an Nth root of 2.
    pub fn with_growth(mut self, growth: f64) -> Self {
        (&mut self).set_from_growth(growth);
        self
    }

    /// Set bins_per_doubling, error, and growth to be consistent with the given value. 
    /// Memory usage is proportional to this value, so this 
    /// is the way to directly control memory usage. 
    /// A bins_per_doubling value of 35 yields 1% error. 
    pub fn with_bins_per_doubling(mut self, bins: usize) -> Self {
        (&mut self).set_from_bins(bins);
        self
    }

    /// Set bins_per_doubling, error, and growth to be consistent with the given error value. 
    /// This is the way to directly control error rate. 
    /// A 1% error corresponds to a bins_per_doubling value of 35. 
    fn set_from_error(&mut self, error: f64) {
        let adjusted_error = 
            if error < 0.00035 { 0.00035 }
            else if error > 0.17 { 0.17 }
            else { error };
        let growth = (1.0 + adjusted_error) / (1.0 - adjusted_error);
        let fractional_bins = 1.0 / growth.log2();
        // Must make the number of bins an integer, then correct the other values.
        self.bins_per_doubling = fractional_bins.ceil() as usize;
        self.growth = 2.0_f64.powf(1.0/(self.bins_per_doubling as f64));
        self.error = (self.growth - 1.0) / (self.growth + 1.0);
    }

    fn set_from_growth(&mut self, growth: f64) {
        let adjusted_growth =
            if growth < 1.00069 { 1.00069 }
            else if growth > 1.4143 { 1.4143 }
            else { growth };
        let fractional_bins = 1.0 / adjusted_growth.log2();
        // Must make the number of bins an integer, then correct the other values.
        self.bins_per_doubling = fractional_bins.ceil() as usize;
        self.growth = 2.0_f64.powf(1.0/(self.bins_per_doubling as f64));
        self.error = (self.growth - 1.0) / (self.growth + 1.0);
    }

    /// Set bins_per_doubling, error, and growth to be consistent with the given bins value. 
    /// This is the way to directly control memory usage. 
    /// The maximum possible memory usage is linearly proportional to bins. 
    /// Actual data usually has much redundancy, autocorrelation, and may obey Zipf's law, 
    /// which will cause much less than the maximum storage in practice. 
    /// 
    /// A bins_per_doubling value of 35 corresponds to 1% error. 
    /// 
    /// Bins will be constrained to be in the range [2,1000].
    fn set_from_bins(&mut self, bins: usize) {
        let adjusted_bins =
            if bins < 2 { 2 }
            else if bins > 1000 { 1000 }
            else { bins };
        self.bins_per_doubling = adjusted_bins;
        self.growth = 2.0_f64.powf(1.0/(self.bins_per_doubling as f64));
        self.error = (self.growth - 1.0) / (self.growth + 1.0);
    }
}

#[derive(Copy, Clone, Debug)]
struct BinInfo {
    pub bin_key : isize,
    pub bin_low_value : f64,
    pub bin_high_value : f64,
    pub bin_weight : f64
}

impl BinInfo {
    pub fn bin_width(&self) -> f64 {
        self.bin_high_value - self.bin_low_value
    }

    /// Fetch the weight of the indicated bin, or zero if it does not exist.
    pub fn get_key_weight(map: &SkipMap<isize, HistogramBin>, key: isize) -> f64 {
        match map.get(&key) {
            Some(bucket) => bucket.weight,
            None => 0.0
        }
    }

    /// Estimate the mode of grouped data using the weights of
    /// the bins immediately below and above.
    /// If those bins are empty, it just takes the midpoint of this bin. 
    /// 
    /// ```text
    ///               /   f1 - f0     \
    ///   Mode = L + |-----------------| h
    ///               \ 2f1 - f0 - f2 /
    /// 
    ///   where:
    ///      L is the lower limit of the modal class
    ///      h is the size of the class interval
    ///      f1 is the frequency of the modal class
    ///      f0 is the frequency of the class preceding the modal class
    ///      f2 is the frequency of the class succeeding the modal class
    /// ```
    pub fn mode_of_grouped_data(&self, bins: &SkipMap<isize, HistogramBin>) -> f64 {
        let f1 = self.bin_weight;
        let f0 = Self::get_key_weight(bins, self.bin_key - 1);
        let f2 = Self::get_key_weight(bins, self.bin_key + 1);
        let h = self.bin_width();
        let l = self.bin_low_value;
        let mode = l + h*(f1 - f0)/(2.0 * f1 - f0 - f2);
        mode
    }
}

/// Maintain a cache of information from which one can estimate
/// the current values for Mode.
struct ModeCache {
    modal_classes : Vec<BinInfo>,
    weight : f64
}

impl ModeCache {
    pub fn new() -> Self {
        ModeCache {
            modal_classes: Vec::new(),
            weight : 0.0
        }
    }

    /// Attempt to update the ModeCache with a new candidate for Mode. 
    /// Return false if no change was made to the mode estimate. 
    pub fn update(&mut self, bin_key : isize,
        bin_low_value : f64,
        bin_high_value : f64,
        bin_weight : f64) -> bool {
        let bin = BinInfo { 
            bin_key : bin_key,
            bin_low_value : bin_low_value,
            bin_high_value : bin_high_value,
            bin_weight : bin_weight
        };
        if self.weight < bin_weight {
            // Replace old mode with a new mode.
            self.modal_classes = vec![bin];
            self.weight = bin_weight;
            true
        }
        else if self.weight == bin_weight {
            // A multi-modal collection of samples.
            self.modal_classes.push(bin);
            true
        }
        else {
            false
        }
    }

    /// Return a list of all modes among the samples. 
    /// If no samples have been added, the list is empty.
    pub fn mode(&self, bins: &SkipMap<isize, HistogramBin>) -> Vec<f64> {
        let mut modes = Vec::new();
        for bin in self.modal_classes.iter() {
            let mode_estimate = bin.mode_of_grouped_data(bins);
            modes.push(mode_estimate);
        }
        modes
    }
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use std::time::{Instant};

    fn assert_close(x: f64, y: f64, epsilon: f64) {
        let delta = (x - y).abs();
        assert!(delta <= epsilon);
    }

    /// Test data with a gap.
    fn gapped_quantogram() -> Quantogram {
        let mut q = Quantogram::new();
        let data = vec![1.0,2.0,3.0,4.0,5.0,95.0,96.0,97.0,98.0,99.0];
        q.add_unweighted_samples(data.iter());
        q
    }

    fn randomish(bottom: f64, top:f64, count: usize) -> Vec<f64> {
        let mut data : Vec<f64> = Vec::new();
        let range = top - bottom;
        let phi: f64 = (1.0 + 5.0_f64.sqrt()) / 2.0;
        for i in 0..count {
            // Pseudo random Weyl sequence based on Golden ratio.
            let pseudo_rand = ((i as f64) * phi) % 1.0_f64;
            // By squaring pseudo_rand we get a non-uniform distribution,
            // which is more interesting.
            let sample = bottom + pseudo_rand * pseudo_rand * range;
            data.push(sample);
        }
        data
    }

    fn assert_median(data: &mut Vec<f64>, accuracy: f64) {
        let mut q = Quantogram::new();
        q.add_unweighted_samples(data.iter());
        let actual_median = q.median().unwrap();
        data.sort_by(|a, b| a.partial_cmp(b).unwrap());
        
        let expected_median = data[data.len() / 2]; // Ignore adjustment for data with an even number of elements.
        let abs_rel_error = (expected_median - actual_median) / actual_median;
        assert!(abs_rel_error <= accuracy);
    }

    /// Add all data and then query a single median.
    fn profile_median(data: &mut Vec<f64>) -> (u128,u128) {
        let t1 = Instant::now();
        let mut q = Quantogram::new();
        q.add_unweighted_samples(data.iter());
        let _actual_median = q.median().unwrap();
        let quantogram_time = t1.elapsed().as_micros();

        let t2 = Instant::now();
        let mut data_copy = Vec::new();
        for sample in data {
          data_copy.push(*sample);
        }
        data_copy.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let _expected_median = data_copy[data_copy.len() / 2]; // Ignore adjustment for data with an even number of elements.
        let naive_time = t2.elapsed().as_micros();
      
        (quantogram_time, naive_time)
    }

    /// Add all data and perform a median after each insertion.
    fn profile_many_medians(data: &mut Vec<f64>) -> (u128,u128) {
        let t1 = Instant::now();
        let mut q = Quantogram::new();
        for sample in data.iter() {
            q.add(*sample);
            let _actual_median = q.median().unwrap();
        }
        let quantogram_time = t1.elapsed().as_micros();

        let t2 = Instant::now();
        let mut data_copy = Vec::new();
        for sample in data {
            data_copy.push(*sample);
            data_copy.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let _expected_median = data_copy[data_copy.len() / 2]; // Ignore adjustment for data with an even number of elements.
        }
        let naive_time = t2.elapsed().as_micros();
      
        (quantogram_time, naive_time)
    }

    fn absolute_relative_error(expected: f64, actual: f64) -> f64 {
        (expected - actual).abs() / expected
    }

    #[test]
    fn test_min() { assert_eq!(gapped_quantogram().min().unwrap(), 1.0); }

    #[test]
    fn test_max() { assert_eq!(gapped_quantogram().max().unwrap(), 99.0); }

    #[test]
    fn test_mean() { assert_eq!(gapped_quantogram().mean().unwrap(), 50.0); }

    #[test]
    fn test_mode_unimodal() { 
        let mut q = Quantogram::new();
        let data = vec![1.0,2.0,3.0,3.0,4.0,4.0,4.0,5.0,6.0];
        q.add_unweighted_samples(data.iter());
        assert_eq!(q.mode(), vec![4.0]); 
    }

    #[test]
    fn test_mode_multimodal() { 
        let mut q = Quantogram::new();
        let data = vec![1.0,2.0,3.0,3.0,3.0,4.0,4.0,4.0,5.0,6.0];
        q.add_unweighted_samples(data.iter());
        assert_eq!(q.mode(), vec![3.0,4.0]); 
    }

    /// Verify that a false outlier mode at zero is accepted by mode but rejected by hsm.
    #[test]
    fn test_hsm() { 
        let mut q = Quantogram::new();
        let data = vec![
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 
            1.0,
            2.0,  2.0,
            3.0,
            4.0,  4.0,  4.0,
            5.0,  5.0,
            6.0,  6.0,  6.0,  6.0,
            7.0,  7.0,  7.0,
            8.0,  8.0,  8.0,  8.0,
            9.0,  9.0,  9.0,  9.0,  9.0,
            10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
            11.0, 11.0, 11.0,
            12.0, 12.0, 12.0,
            13.0,
            14.0
        ];
        q.add_unweighted_samples(data.iter());
        assert_eq!(q.hsm(), Some(10.0)); 
        assert_eq!(q.mode(), vec![0.0]); 
    }

    #[test]
    fn test_count() { 
        assert_eq!(Quantogram::new().count(), 0); 
        assert_eq!(gapped_quantogram().count(), 10); 
    }

    #[test]
    fn test_median() { assert_eq!(gapped_quantogram().median().unwrap(), 5.0); }


    #[test]
    fn test_unweighted_standard_deviation() { 
        let mut q = Quantogram::new();

        q.add_unweighted_samples(vec![
            85, 86, 100, 76, 81, 93, 84, 99, 71, 69, 93, 85, 81, 87, 89
        ].iter());

        let expected_stddev = 8.698403429493382;
        assert_close(q.stddev().unwrap(), expected_stddev, 0.00000000001);
    }

    #[test]
    fn test_weighted_standard_deviation() { 
        let mut q = Quantogram::new();

        q.add_unweighted_samples(vec![
            86, 100, 76, 93, 84, 99, 71, 69, 93, 87, 89
        ].iter());
        q.add_weighted(85.0, 2.0);
        q.add_weighted(81.0, 2.0);

        let expected_stddev = 8.69840342949338;
        assert_close(q.stddev().unwrap(), expected_stddev, 0.00000000001);
    }

    /// Removing values may mess up the standard deviation, but this case passes.
    #[test]
    fn test_unweighted_standard_deviation_with_remove() { 
        // Include spurious values of 123 and 150, then remove them. 
        let mut q = Quantogram::new();
        q.add_unweighted_samples(vec![
            123, 150, 85, 86, 100, 76, 81, 93, 84, 99, 71, 69, 93, 85, 81, 87, 89
        ].iter());
        q.remove(150.0);
        q.remove(123.0);

        let expected_stddev = 8.698403429493382;
        assert_close(q.stddev().unwrap(), expected_stddev, 0.00000000001);
    }

    #[test]
    fn test_quantile() { 
        assert_eq!(gapped_quantogram().quantile(0.0).unwrap(), 1.0); 
        assert_eq!(gapped_quantogram().quantile(0.75).unwrap(), 96.0); 
        assert_eq!(gapped_quantogram().quantile(1.0).unwrap(), 99.0); 
    }

    #[test]
    fn test_quantile_at() { 
        let q_mean = |r: Option<(f64,f64)>| { (r.unwrap().0 + r.unwrap().1) / 2.0 };
        let mut q = Quantogram::new();
        let data: Vec<f64> = vec![0,1,2,3,4,5,6,7,8,9,10]
            .into_iter()
            .map(|x| x as f64)
            .collect();
        q.add_unweighted_samples(data.iter());
        assert_eq!(0.5, q_mean(q.quantile_at(5.0)));
        assert_eq!(17.0/22.0, q_mean(q.quantile_at(8.0)));
    }

    #[test]
    fn test_fussy_median() { 
        assert_eq!(gapped_quantogram().fussy_quantile(0.5, 2.0).unwrap(), 50.0); 
    }

    #[test]
    fn test_remove() { 
        let mut q = gapped_quantogram();
        q.remove(99.0);
        q.remove(98.0);
        assert_eq!(q.median().unwrap(), 4.0); 
    }

    #[test]
    fn test_quantiles_with_negatives() { 
        let mut q = Quantogram::new();
        let data = vec![-1.0,-2.0,-3.0,-4.0,-5.0,-95.0,-96.0,-97.0,-98.0,-99.0, 0.0, 1.0,2.0,3.0,4.0,5.0,95.0,96.0,97.0,98.0,99.0];
        q.add_unweighted_samples(data.iter());

        assert_eq!(q.median().unwrap(), 0.0); 
        assert_eq!(q.quantile(0.4).unwrap(), -2.0); 
        assert_eq!(q.quantile(0.6).unwrap(), 2.0); 
    }

    #[test]
    fn test_weighted_quantile() { 
        let mut q = Quantogram::new();
        let data = vec![1.0,2.0,3.0,4.0,5.0,95.0,96.0,97.0,98.0,99.0];
        q.add_unweighted_samples(data.iter());
        q.add_weighted(0.0, 6.0);

        assert_eq!(q.median().unwrap(), 2.0); 
    }

    #[test]
    fn test_coeff_of_range() { 
        let mut q = Quantogram::new();
        q.add_unweighted_samples(vec![
            123, 150, 85, 86, 100, 76, 81, 93, 84, 99, 71, 69, 93, 85, 81, 87, 89
        ].iter());

        let expected_coeff = (150.0 - 69.0)/(150.0 + 69.0);
        assert_close(q.coeff_of_range().unwrap(), expected_coeff, 0.00000000001);
    }

    #[test]
    fn test_coeff_of_quartile_deviation() { 
        let mut q = Quantogram::new();
        q.add_unweighted_samples(vec![
            1,2,3,4,4,5,6,7,7,8,9,10,10,11,12,13
        ].iter());

        let expected_coeff = (10.0 - 4.0)/(10.0 + 4.0);
        assert_close(q.coeff_of_quartile_dev().unwrap(), expected_coeff, 0.00000000001);
    }

    #[test]
    fn test_quartiles() { 
        let mut q = Quantogram::new();
        q.add_unweighted_samples(vec![
            1,2,3,4,4,5,6,7,7,8,9,10,10,11,12,13
        ].iter());

        assert_close(q.q1().unwrap(), 4.0, 0.00000000001);
        assert_close(q.q3().unwrap(), 10.0, 0.00000000001);
    }

    #[test]
    fn test_coeff_of_stddev() { 
        let mut q = Quantogram::new();
        q.add_unweighted_samples(vec![
            85, 86, 100, 76, 81, 93, 84, 99, 71, 69, 93, 85, 81, 87, 89
        ].iter());

        let expected_stddev =  8.698403429493382;
        let expected_mean   = 85.266666666666667;
        assert_close(
            q.coeff_of_stddev().unwrap(), 
            expected_stddev / expected_mean, 
            0.00000000001
        );
        // Also test the coeff_of_variation
        assert_close(
            q.coeff_of_variation().unwrap(), 
            100.0 * expected_stddev / expected_mean, 
            0.00000000001
        );
    }

    #[test]
    fn test_finite() { 
        let mut q = Quantogram::new();
        q.add(f64::NAN);
        assert_eq!(q.finite(), 0.0);

        let data = vec![1.0,2.0,3.0,4.0];
        q.add_unweighted_samples(data.iter());
        assert_eq!(q.finite(), 0.8); 
    }

    #[test]
    fn test_nan() { 
        let mut q = Quantogram::new();
        q.add(f64::NAN);
        assert_eq!(q.nan(), 1.0);

        let data = vec![1.0,2.0,3.0,4.0];
        q.add_unweighted_samples(data.iter());
        assert_eq!(q.nan(), 0.2); 
    }

    #[test]
    fn test_size() { 
        let mut q = Quantogram::new();

        // Despite sparseness, storage can go no lower than 7, because we always record counts for NANs, zeroes, etc.
        assert_eq!(q.size(), 7); 
        
        // Adding the first value should add a bin to a coarse and a fine SkipMap.
        q.add(1.0);
        assert_eq!(q.size(), 9); 

        // Adding the same value should not increase the storage used. 
        q.add(1.0);
        assert_eq!(q.size(), 9); 

        // Adding a new value that is in the same coarse bin but a different fine bin should increase the storage by one, not two. 
        q.add(1.5);
        assert_eq!(q.size(), 10); 
    }

    /// Verify that storage grows with the logarithm of the number of samples. 
    #[test]
    fn test_sparsity() { 
        let mut q = Quantogram::new();
        // A million is about 20 powers of 2. 
        // Default histogram has 35 bins per power of two.
        // This means that we should have 7 default bins, 20 coarse bins and at most 700 fine bins. 
        for i in 1..1000000 {
            q.add(i as f64);
        }
        // Actual size is 577.
        assert!(q.size() < 600);

        assert!(absolute_relative_error(q.median().unwrap(), 500000.0) <= 0.01);
    }

    #[test]
    fn test_with_constrained_range() { 
        let mut q = QuantogramBuilder::new()
            .with_smallest_power(-10)
            .with_largest_power(10)
            .build();
        let data = vec![
            0.0001, 0.00056, // will underflow and become zeroes
            0.002, 16.0, 1000.0, 
            6543.0, 15000.0, 2400000.0 // will overflow and become Infinity
        ];
        q.add_unweighted_samples(data.iter());
        assert_eq!(q.min().unwrap(), 0.0); 
        assert_eq!(q.finite(), 0.625); 
        assert_eq!(q.zero(), 0.25); 
        assert_eq!(q.median().unwrap(), 0.002); 
    }

    #[test]
    fn test_accuracy() {
        let mut data = randomish(-100.0, 2000000.0, 100000);
        assert_median(&mut data, 0.01);
    }

    #[test]
    fn test_insert_heavy_speed() {
        // Adding more data to a naive quantile should take longer compared to using a histogram
        // because the naive approach has an N log N sort. 
        // This test verifies that ther is an improvement in relative performance
        // between the two as the number of samples increases.

        let mut short_data = randomish(-10000.0, 10000.0, 10000);
        let mut long_data = randomish(-10000.0, 10000.0, 1000000);
        let (short_quant_time, short_naive_time) = profile_median(&mut short_data);
        let (long_quant_time, long_naive_time) = profile_median(&mut long_data);
        let short_ratio = (short_quant_time as f64) / (short_naive_time as f64);
        let long_ratio = (long_quant_time as f64) / (long_naive_time as f64);
        println!("Short Ratio: {}  Long ratio: {}", short_ratio, long_ratio);
        // assert!(false);
        assert!(short_ratio > 1.4 * long_ratio);
    }

    #[test]
    fn test_median_heavy_speed() {
        let mut short_data = randomish(-10000.0, 10000.0, 2000);
        let mut long_data = randomish(-10000.0, 10000.0, 20000);
        let (short_quant_time, short_naive_time) = profile_many_medians(&mut short_data);
        let (long_quant_time, long_naive_time) = profile_many_medians(&mut long_data);
        let short_ratio = (short_quant_time as f64) / (short_naive_time as f64);
        let long_ratio = (long_quant_time as f64) / (long_naive_time as f64);
        println!("Short Ratio: {}  Long ratio: {}", short_ratio, long_ratio);
        // assert!(false);
        assert!(short_ratio > 8.0 * long_ratio);
    }

    #[test]
    fn test_basics() {
        let mut q = Quantogram::new();
        q.add(10.0);
        q.add(40.0);
        q.add(20.0);
        q.add(30.0);
        q.add(50.0);
    
        assert_eq!(q.count(), 5);
        assert_eq!(q.min().unwrap(), 10.0);
        assert_eq!(q.max().unwrap(), 50.0);
        assert_eq!(q.mean().unwrap(), 30.0);
        assert_eq!(q.median().unwrap(), 30.0);
        assert_eq!(q.quantile(0.75).unwrap(), 40.0);

        q.remove(10.0);
        q.remove(20.0);
        assert_eq!(q.mean().unwrap(), 40.0);
    }

    #[test]
    fn quantile_at_nan() {
        let mut q = Quantogram::new();
        q.add(10.0);
        assert!(q.quantile_at(f64::NAN).is_none()); 
    }

}

