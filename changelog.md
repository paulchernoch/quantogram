# Quantogram Release Notes


**Release 0.4.2**:

Remove debug logging statement from hsm.

**Release 0.4.1**:

Fix half-sample mode (hsm) edge case where there is no true mode.

**Release 0.4.0**:

Added these basic statistics:

 - **q1** (1st quartile)
 - **q3** (3rd quartile)

Added these measures of dispersion: 

 - **range** (spread between min and max)
 - **iqr** (interquartile range: Q3 - Q1)        
 - **quartile_deviation** (semi-interquartile range, half of iqr))
 - **coeff_of_range** ((max - min) / (max + min))
 - **coeff_of_quartile_dev** ((Q3-Q1)/(Q3+Q1))
 - **coeff_of_stddev** (standard deviation divided by mean)
 - **coeff_of_variation**

**Release 0.3.0**:

 - Added the **half-sample mode** measure, `Quantogram::hsm`, an estimate of the mode that is resistant to outliers. Based on the Robertson-Cryer half-sample mode.
 - If adding the same sample multiple times and a bin has only one distinct value that matches that accurate value, keep using the accurate value and do not switch to using the bin's midpoint value.

  
**Release 0.2.0**:

 - Added Quantogram::variance
 - Added Quantogram::stddev (standard deviation)
 - Quantogram::add_unweighted_samples now accepts iterators over any value that can convert to an f64 via the Into trait.