# LAMBDA 4.0

The runtime comparisons call the externally installed LAMBDA 4.0 toolbox.
Obtain the toolbox and its terms from the authors at Delft University of Technology.

[LAMBDA download page](https://www.tudelft.nl/en/ceg/about-faculty/departments/geoscience-remote-sensing/research/lambda/lambda/)

The counted reference routine `matlab/metrics/lambda4_reduce_covariance_counted.m`
follows the algorithms in LAMBDA 4.0 `decomposeLtDL`, `transformZ`, and
`computeIGT_row`, with explicit arithmetic and event counters. The GPS covariance
generator extends the model in LAMBDA 4.0 `RUN_example_3` to three frequencies.
These sources are credited to Lotfi Massarweh and the Geoscience & Remote Sensing
department at TU Delft, whose toolbox headers carry a June 2024 copyright notice.
The MIT license covers the original contributions in this repository. The external
toolbox remains subject to its own terms and copyright notices.

Massarweh, L., Verhagen, S., and Teunissen, P. J. G. (2025).
New LAMBDA toolbox for mixed-integer models: estimation and evaluation.
GPS Solutions, 29, 14. https://doi.org/10.1007/s10291-024-01738-z
