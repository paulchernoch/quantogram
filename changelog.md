# Quantogram Release Notes

**Release 0.3.0**:

 - Added the **half-sample mode** measure, `Quantogram::hsm`, an estimate of the mode that is resistant to outliers. Based on the Robertson-Cryer half-sample mode.
 - If adding the same sample multiple times and a bin has only one distinct value that matches that accurate value, keep using the accurate value and do not switch to using the bin's midpoint value.

  
**Release 0.2.0**:

 - Added Quantogram::variance
 - Added Quantogram::stddev (standard deviation)
 - Quantogram::add_unweighted_samples now accepts iterators over any value that can convert to an f64 via the Into trait.