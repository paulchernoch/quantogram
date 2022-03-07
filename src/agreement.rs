

/// Records how much in agreement are a measurement and a proxy for that measurement.
/// If an inexpensive proxy is close enough to an expensive-to-compute measure then
/// the proxy may be used in its place.
#[derive(Clone,Copy,Debug)]
pub struct Agreement {
    /// Cumulative counts of how similar a measure is with its proxy.
    ///   - Index 0 is full count of comparisons made. 
    ///   - Index 1 is count of measurements where abs(proxy - measure) / measure <= 1/2
    ///   - Index 2 is count of measurements where abs(proxy - measure) / measure <= 1/4
    ///   - Index 3 is count of measurements where abs(proxy - measure) / measure <= 1/8
    ///   - All the way up to index 30, for <= 1/2^30
    ///   - Index 31 is for count of when proxy = measure exactly 
    cumulative_similarity: [usize;32]
}

impl Agreement {
    pub fn new() -> Self {
        Agreement {
            cumulative_similarity: [0;32]
        }
    }

    /// Compare measure to proxy and add to histogram.
    pub fn add(&mut self, measure: f64, proxy: f64) {
        self.cumulative_similarity[0] += 1;
        if measure == 0.0 {
            if proxy == 0.0 {
                for i in 1..32 {
                    self.cumulative_similarity[i] += 1;
                }
            }
        }
        else {
            let abs_rel_error = (proxy - measure).abs() / measure;
            let mut curr_error = 0.5;
            for i in 1..31 {
                if abs_rel_error <= curr_error {
                    self.cumulative_similarity[i] += 1;
                }
                else {
                    break;
                }
                curr_error /= 2.0;
            }
            if abs_rel_error == 0.0 {
                self.cumulative_similarity[31] += 1;
            }
        }
    }

    /// How often does the proxy match the measure exactly, as a fraction in [0,1].
    pub fn matches_exactly(&self) -> f64 {
        if self.cumulative_similarity[0] == 0 { 0.0 }
        else {
            self.cumulative_similarity[31] as f64 / self.cumulative_similarity[0] as f64
        }
    }

    /// How often does the proxy match the measure with error <= 1/2^power ?
    pub fn matches_approximately(&self, power: usize) -> f64 {
        if self.cumulative_similarity[0] == 0 { 0.0 }
        else if power > 30 {
            self.cumulative_similarity[31] as f64 / self.cumulative_similarity[0] as f64
        }
        else {
            self.cumulative_similarity[power] as f64 / self.cumulative_similarity[0] as f64
        }
    }

    /// Get the highest power such that the error is 1/2^power or less with 
    /// the given likelihood or better.
    /// If exact matches are that likely, return 31.
    pub fn get_matching_power(&self, likelihood: f64) -> usize {
        if self.cumulative_similarity[0] == 0 { 0 }
        else if self.matches_exactly() >= likelihood { 31 }
        else {
            let denominator = self.cumulative_similarity[0] as f64;
            // TODO: Could use binary search here.
            for power in 1..31 {
                let ratio = self.cumulative_similarity[power] as f64 / denominator;
                if likelihood < ratio { 
                    return power - 1;
                }    
            }
            30          
        } 
    }

    /// Decides if the proxy agrees with the measure, having an absolute relative error <= 1/2^power
    /// at least as often as likelihood (a fraction in the range [0,1]).
    pub fn agrees(&self, power: usize, likelihood: f64) -> bool {
        self.matches_approximately(power) >= likelihood
    }
}

