use std::cmp::Ordering;

// Contains: 
// 
//   half_sample_mode - Computes the half-sample mode, an estimate of the mode for weighted data.
//   HalfSampleInterval - internal struct
//   HalfSampleModeCache - Caches prior result of half_sample_mode to focus recalculation of mode on a subset of data, a time optimization. 

/// Internal structure for Half Sample Mode calculation.
#[derive(Copy, Clone, Debug)]
struct HalfSampleInterval {
    /// Index of the first bin in rthe interval
    pub start_index: usize,

    /// Sample value for first bin.
    pub start_sample: f64,

    /// Count of the bins in the interval
    pub bins: usize,

    /// Sum of weights from all bins in the interval
    pub cume_weight: f64,

    /// Width of interval is the difference between the lowest bin sample value 
    /// and the highest bin sample value in the interval.
    pub width: f64
}

impl HalfSampleInterval {
    pub fn new(index: usize, sample: f64, weight: f64) -> Self {
        HalfSampleInterval {
            start_index: index,
            start_sample: sample,
            bins: 1,
            cume_weight: weight,
            width: 0.0
        }
    }

    pub fn reset(&mut self, index: usize, sample: f64, weight: f64) {
        self.start_index = index;
        self.start_sample = sample;
        self.bins = 1;
        self.cume_weight = weight;
        self.width = 0.0;
    }

    pub fn empty() -> Self {
        Self::new(0, 0.0, 0.0)
    }

    pub fn is_empty(&self) -> bool {
        self.cume_weight == 0.0
    }

    pub fn add(&mut self, sample: f64, weight: f64) {
        self.bins += 1;
        self.cume_weight += weight;
        self.width = sample - self.start_sample;
    }

    pub fn is_underweight(&self, minimum_weight: f64) -> bool {
        self.cume_weight < minimum_weight
    }

    /// Check if this interval is denser than the other interval. 
    pub fn is_denser_than(&self, other: &Self, minimum_weight: f64) -> bool {
        if other.is_empty() && !self.is_empty() { true }
        else if self.is_underweight(minimum_weight) { false }
        else if other.is_underweight(minimum_weight) { true }
        else if self.width < other.width { true }
        else if self.width > other.width { false }
        else { self.cume_weight > other.cume_weight }
    }

}

/// Estimate the mode using the half-sample mode algorithm.
/// 
/// More resistant than the mode to outliers and contamination (noise). 
/// 
/// Based on the mode estimator by Robertson and Cryer (1974), and an 
/// algorithm described in "On a Fast, Robust Estimator of the Mode" by David Bickel and Rudolf Fruhwirth.
/// That algorithm is applied to raw samples, whereas this is applied 
/// to already histogrammed values that are weighted.
/// 
/// NOTE: This has worse than quadratic performance, possibly N^2 Log N.
/// 
/// # Arguments
/// 
/// * `samples` - List of sample values, sorted in ascending order.
/// * `weights` - List of positive weights, corresponding to samples.
/// * `total_weight` - Sum of the weights.
pub fn half_sample_mode(samples: &[f64], weights: &[f64], total_weight: f64) -> Option<f64> {
    if total_weight == 0.0 { None }
    else {
        let count = samples.len();
        if count < 5 {
            let mut max_weight = 0.0;
            let mut samples_with_max_weight: Vec<f64> = Vec::new();
            for (&sample, &weight) in samples.iter().zip(weights) {
                if weight > max_weight {
                    max_weight = weight;
                    samples_with_max_weight = vec![sample];
                }
                else if weight == max_weight {
                    samples_with_max_weight.push(sample);
                }
            }
            return Some(samples_with_max_weight.iter().sum::<f64>() / samples_with_max_weight.len() as f64);
        } 
        let mut cume_target = total_weight / 2.0;
        let mut start_index = 0;
        let mut limit_index = samples.len();
        let mut densest = HalfSampleInterval::empty();
        while (limit_index - start_index) >= 2 { 
            densest = HalfSampleInterval::empty();
            for interval_start in start_index..limit_index {
                let mut current_interval = HalfSampleInterval::empty();
                let mut made_weight = false;
                for index in interval_start..limit_index {
                    if current_interval.is_empty() {
                        current_interval.reset(index, samples[index], weights[index]);
                    }
                    else {
                        current_interval.add(samples[index], weights[index]);
                    }
                    if current_interval.is_denser_than(&densest, cume_target) {
                        densest = current_interval;
                    }
                    if !current_interval.is_underweight(cume_target) { 
                        made_weight = true;
                        break; 
                    }     
                }
                // Even going all the way to limit_index did not accumulate enough
                // weight to equal cume_target, so we can be done with the middle loop.
                if !made_weight {
                    break;
                }
            }
            start_index = densest.start_index;
            limit_index = start_index + densest.bins;
            cume_target = cume_target / 2.0;
        }
        match densest.bins {
            0 => None,
            1 => Some(samples[densest.start_index]),
            2 => {
                // TODO: Instead of taking the weighted average of two points, should we take the maximum?
                let i = densest.start_index;
                let weight_sum = weights[i] + weights[i+1];
                if weight_sum <= 0.0 {
                    None
                }
                else {
                    let weighted_mean = (samples[i] * weights[i] + samples[i+1] * weights[i+1]) / weight_sum;
                    Some(weighted_mean)
                }
            },
            3 => {
                // TODO: Instead of taking the weighted average of three points, should we take the maximum?
                let i = densest.start_index;
                let weight_sum = weights[i] + weights[i+1] + weights[i+2];
                if weight_sum <= 0.0 {
                    None
                }
                else {
                    let weighted_mean = (
                        samples[i] * weights[i] 
                      + samples[i+1] * weights[i+1]
                      + samples[i+2] * weights[i+2]) 
                      / weight_sum;
                    Some(weighted_mean)
                }
            },
            _ => None
        }
    }
}

/// Permit the caching and reuse of mode values computed using half-sample mode
/// in order to enhance performance. 
/// 
/// half_sample_mode (hsm) is a costly function, worse than quadratic.
/// 
/// Cache the last computed value of the half-sample mode and 
/// use it to restrict the range of samples to be reazalyzed
/// when you compute a new mode after one or more samples have been inserted. 
/// 
/// If the new sample is in the same neighborhood as the previous mode, 
/// assume that the new mode must be in the same neighborhood and
/// only make computations on samples in that neighborhood.
/// 
/// The cache is invalidated if: 
///   - the new sample is out of the neighborhood
///   - the total_weight has changed too much
///   - the number of bins has changed too much
#[derive(Copy, Clone, Debug)]
pub struct HalfSampleModeCache {
    /// Previous mode computed by hsm, or none if invalidated.
    mode: Option<f64>,

    /// Most recent sample added to date.
    latest_sample: f64,

    /// Number of bins present during previous full recalculation of mode.
    bins: usize,

    /// Total weight of all samples during previous full recalculation of mode.
    total_weight: f64,

    /// Number of consecutive bins considered near to the current mode, half below and half above.
    neighborhood: usize,

    /// If the number of bins has changed by more than this many, invalidate the cache.
    allowed_bin_change: usize,

    /// If the relative change in total_weight exceeds this, invalidate the cache.
    allowed_weight_change_fraction: f64,

    allowed_outside_neighborhood: usize,

    outside_neighborhood_counter: usize,

    cache_hits: usize
}

impl HalfSampleModeCache {
    pub fn new(neighborhood: usize, allowed_bin_change: usize, allowed_weight_change_fraction: f64) -> Self {
        HalfSampleModeCache {
            mode: None,
            latest_sample: f64::NAN,
            bins: 0,
            total_weight: 0.0,
            neighborhood: neighborhood,
            allowed_bin_change: allowed_bin_change,
            allowed_weight_change_fraction: allowed_weight_change_fraction,
            allowed_outside_neighborhood: 12,
            outside_neighborhood_counter: 0,
            cache_hits: 0
        }
    }

    pub fn new_default() -> Self {
        HalfSampleModeCache {
            mode: None,
            latest_sample: f64::NAN,
            bins: 0,
            total_weight: 0.0,
            // A neighborhood of 70 is chosen so that if a growth factor of 1.02 is used
            // in a Quantogram for 1% error, we will range over data in the 
            // range [1/2 x mode, 2 x mode] because there are then 35 bins per doubling,
            // which should capture sufficient range.
            neighborhood: 70,
            allowed_bin_change: 10,
            allowed_weight_change_fraction: 0.01,
            allowed_outside_neighborhood: 12,
            outside_neighborhood_counter: 0,
            cache_hits: 0
        }
    }

    /// Record the latest sample added to data.
    /// This must be called prior to calling hsm so that we know if the new value is near the prior mode.
    pub fn record(&mut self, sample: f64) {
        self.latest_sample = sample;
    }

    pub fn hits(&self) -> usize { self.cache_hits }

    fn full_recalc(&mut self, samples: &[f64], weights: &[f64], total_weight: f64) -> Option<f64> {
        self.bins = samples.len();
        self.total_weight = total_weight;
        self.mode = half_sample_mode(samples, weights, total_weight);
        self.outside_neighborhood_counter = 0;
        self.mode
    }

    /// Cribbed from the nightly feature on f64 to order floats.
    /// TODO: Once that feature becomes stable, should switch to using it.
    /// See #[unstable(feature = "total_cmp", issue = "72599")]
    #[inline]
    pub fn total_cmp(this: &f64, other: &f64) -> Ordering {
        let mut left = this.to_bits() as i64;
        let mut right = other.to_bits() as i64;

        // In case of negatives, flip all the bits except the sign
        // to achieve a similar layout as two's complement integers
        //
        // Why does this work? IEEE 754 floats consist of three fields:
        // Sign bit, exponent and mantissa. The set of exponent and mantissa
        // fields as a whole have the property that their bitwise order is
        // equal to the numeric magnitude where the magnitude is defined.
        // The magnitude is not normally defined on NaN values, but
        // IEEE 754 totalOrder defines the NaN values also to follow the
        // bitwise order. This leads to order explained in the doc comment.
        // However, the representation of magnitude is the same for negative
        // and positive numbers – only the sign bit is different.
        // To easily compare the floats as signed integers, we need to
        // flip the exponent and mantissa bits in case of negative numbers.
        // We effectively convert the numbers to "two's complement" form.
        //
        // To do the flipping, we construct a mask and XOR against it.
        // We branchlessly calculate an "all-ones except for the sign bit"
        // mask from negative-signed values: right shifting sign-extends
        // the integer, so we "fill" the mask with sign bits, and then
        // convert to unsigned to push one more zero bit.
        // On positive values, the mask is all zeros, so it's a no-op.
        left ^= (((left >> 63) as u64) >> 1) as i64;
        right ^= (((right >> 63) as u64) >> 1) as i64;

        left.cmp(&right)
    }

    /// Recalculate the half-sample mode but only examine a neighborhood around the previous
    /// estimate for mode, not the whole data set.
    fn partial_recalc(&mut self, prev_mode_position: usize, samples: &[f64], weights: &[f64], total_weight: f64) -> Option<f64> {
        self.cache_hits += 1;
        self.bins = samples.len();
        self.total_weight = total_weight;
        let half_neighborhood = self.neighborhood / 2;

        let start_index = if prev_mode_position < half_neighborhood { 0 }
        else { prev_mode_position - half_neighborhood };

        let stop_index = if prev_mode_position + half_neighborhood >= samples.len() { samples.len() }
        else { prev_mode_position + half_neighborhood };

        self.mode = half_sample_mode(&samples[start_index..stop_index], &weights[start_index..stop_index], total_weight);
        self.mode
    }

    /// Compute the half-sample mode of weighted data. 
    /// 
    /// A full recalculation or a partial recalculation will be performed, depending on
    /// whether the potential for the mode to have drifted far seems likely or not.
    /// 
    /// Before invoking this, you must call record(sample) to set the latest sample value.
    /// This is so that you can add several samples in between calls to hsm if necessary, then call hsm.
    /// The adder of samples should not need to know when or if hsm is to be called soon.
    pub fn hsm(&mut self, samples: &[f64], weights: &[f64], total_weight: f64) -> Option<f64> {
        let new_bins = samples.len();
        let do_full_recalc = self.mode.is_none()
          || new_bins < 64 // For fewer samples, no need to cut corners because fast enough
          || new_bins - self.bins > self.allowed_bin_change
          || (total_weight - self.total_weight) / self.total_weight > self.allowed_weight_change_fraction
        ;
        let prev_mode = self.mode.unwrap_or_default();
        if !do_full_recalc && prev_mode.is_finite() {
            let prev_mode_position = match samples.binary_search_by(|probe| Self::total_cmp(probe, &prev_mode)) {
                Result::Ok(pos) => pos,
                Err(pos) => pos
            };
            let newest_sample_position = match samples.binary_search_by(|probe| Self::total_cmp(probe, &self.latest_sample)) {
                Result::Ok(pos) => pos,
                Err(pos) => pos
            };
            // If the added sample is too far from the center of the neighborhood, we need a full recalc.
            let spread = usize::max(prev_mode_position, newest_sample_position) - usize::min(prev_mode_position, newest_sample_position);
            let third_neighborhood = self.neighborhood / 3;
            let half_neighborhood = self.neighborhood / 2;
            if spread > half_neighborhood && self.outside_neighborhood_counter < self.allowed_outside_neighborhood {
                // If the latest sample is outside the neighborhood, is unlikely to move the mode. 
                // Increment a counter so that we are forced to invalidate the cache 
                // if too many such samples are received. 
                self.outside_neighborhood_counter += 1;
                // Reuse previous value with no recalcuation.
                self.cache_hits += 1;
                self.mode
            }
            else if spread > third_neighborhood {
                if self.outside_neighborhood_counter >= self.allowed_outside_neighborhood {             
                    self.full_recalc(samples, weights, total_weight)
                }
                else {
                    self.outside_neighborhood_counter += 1;
                    self.partial_recalc(prev_mode_position, samples, weights, total_weight)
                }
            }
            else {
                self.partial_recalc(prev_mode_position, samples, weights, total_weight)
            }
        }
        else {
            self.full_recalc(samples, weights, total_weight)
        }
    }

    /// Clear the cache, forcing the next call to hsm to do a full recalc.
    pub fn clear(&mut self) {
        self.mode = None;
        self.latest_sample = f64::NAN;
        self.bins = 0;
        self.total_weight = f64::NAN;
        self.cache_hits = 0;
        self.outside_neighborhood_counter = 0;
    }
}


#[cfg(test)]
mod tests {
    use std::time::{Duration, Instant};
    use super::*;

    #[test]
    fn hsm_cache_performance() {
        let no_cache_duration = simulate_convergence(1200, false).as_secs_f64();
        let with_cache_duration = simulate_convergence(1200, true).as_secs_f64();
        let speedup = no_cache_duration / with_cache_duration;
        println!("
        With cache:    {:.2} sec
        Without Cache: {:.2} sec
        Speedup:       {:.3}x
        ", with_cache_duration, no_cache_duration, speedup);
        assert!(speedup > 8.0);
    }

    /// Golden ratio based pseudo-random number sequence.
    ///  ```text
    ///    φ = (1 + √5) / 2
    ///    G(n) = (n · φ) mod 1
    /// ```
    fn golden_ratio_sequence(n: usize) -> f64 {
        let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
        (n as f64 * phi) % 1.0
    }

    fn insert_sorted(value: f64, samples: &mut Vec<f64>) {
        let pos = match samples.binary_search_by(|probe| HalfSampleModeCache::total_cmp(probe, &value)) {
            Result::Ok(pos) => pos,
            Err(pos) => pos
        };
        samples.insert(pos, value);
    }

    /// Simulate a series of estimates of a value that steadily converge. 
    /// That means that a modal peak will grow near the true solution.
    fn simulate_convergence(iterations: usize, use_cache: bool) -> Duration {
        let true_value = 100.0;
        let start_value = 0.0;
        let start_error = start_value - true_value;
        let mut values: Vec<f64> = Vec::new();
        let mut weights: Vec<f64> = Vec::new();
        let mut weight_sum = 0.0;
        let mut cache = HalfSampleModeCache::new_default();
        let start_time = Instant::now();
        for i in 1..=iterations {
            let error = start_error / (i as f64).sqrt();
            let sign = if golden_ratio_sequence(i) < 0.3 { -1.0 } else { 1.0 };
            let next_value = true_value + sign * error;
            // Must ensure values remain sorted.
            insert_sorted(next_value, &mut values);
            // values.push(next_value);
            weights.push(1.0);
            weight_sum += 1.0;
            if use_cache {
                cache.record(next_value);
                let _mode = cache.hsm(&values, &weights, weight_sum);
            }
            else {
                let _mode = half_sample_mode(&values, &weights, weight_sum);
            }
        }
        if use_cache {
            println!("  Cache hits: {}. Cache misses: {}", cache.hits(), iterations - cache.hits());
        }
        start_time.elapsed()
    }

}