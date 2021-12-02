# simples

A rust library for smoothing and simplification of N-dimensional linestrings (as `&[nalgebra::Point]`s).

Intended for cases where the endpoints of the linestrings cannot be moved
(for example, if they touch other objects at those locations).

Currently supports:

- Simplification
  - Resampling at arbitrary distances
  - Ramer-Douglass-Peucker
  - Visvalingam-Whyatt
- Smoothing
  - Moving average
  - Gaussian
  - A Kernel trait for implementing your own kernels to drop in

## To do

- Mapping old points on to resampled points
- Traits for smoothable/ resampleable?
