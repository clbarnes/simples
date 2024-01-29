//! Smooth linestrings.
//!
//! Linestrings are smoothed if they keep the same number of points, but move them around.
use crate::Precision;
use nalgebra::{distance_squared, Point};
use std::cmp::Ordering;
use std::collections::hash_map::Entry;
use std::collections::HashMap;

fn mean<const D: usize>(points: &[Point<Precision, D>]) -> Point<Precision, D> {
    // can't do a straight reduce because of ownership stuff
    let mut iter = points.iter();
    let first = *iter.next().unwrap();
    points.iter().fold(first, |prev, curr| prev + curr.coords) / points.len() as Precision
}

fn weighted_mean<const D: usize>(
    point_weights: &[(Point<Precision, D>, Precision)],
) -> Point<Precision, D> {
    let (point, weight) = point_weights
        .iter()
        .map(|(point, weight)| (point * *weight, *weight))
        .reduce(|prev, next| (prev.0 + next.0.coords, prev.1 + next.1))
        .unwrap();
    point / weight
}

/// Smooth line using a moving average with `2*width + 1` points centred on the point of interest.
/// At the ends, uses a smaller window.
///
/// Should probably only be used on a line already resampled with same-length gaps.
pub fn smooth_moving_average<const D: usize>(
    line: &[Point<Precision, D>],
    width: usize,
) -> Vec<Point<Precision, D>> {
    if line.len() <= 2.max(width) {
        return line.to_vec();
    }
    let mut out = vec![*line.first().unwrap()];
    for this_width in 1..width {
        let this_window = this_width * 2 + 1;
        out.push(mean(&line[0..this_window]));
    }
    for window in line.windows(2 * width + 1) {
        out.push(mean(window));
    }
    for this_width in 1..width {
        let this_window = this_width * 2 + 1;
        out.push(mean(&line[(line.len() - this_window)..]));
    }
    out.push(*line.last().unwrap());

    out
}

/// Structs which can be use as a smoothing kernel.
pub trait Kernel {
    /// If a point is `dist` away from the point of interest, how much should we care about its position?
    ///
    /// None if too low to care about or otherwise invalid.
    /// Weights do not need to sum up to 1.
    fn weigh_dist(&self, dist: Precision) -> Option<Precision>;

    /// If the squared distance from the point of interest to another point is `dist2`, how much should we care about its position?
    ///
    /// None if too low to care about or otherwise invalid.
    /// Weights do not need to sum up to 1.
    fn weigh_dist2(&self, dist2: Precision) -> Option<Precision> {
        self.weigh_dist(dist2.sqrt())
    }

    /// The weight of the point of interest.
    fn at_center(&self) -> Precision;
}

/// Weight points by how far they are from the point of interest in a linear fashion.
#[derive(Copy, Clone, Debug)]
pub struct Linear {
    max_dist: Precision,
}

impl Linear {
    pub fn new(max_dist: Precision) -> Self {
        Self { max_dist }
    }
}

impl Kernel for Linear {
    fn weigh_dist(&self, dist: Precision) -> Option<Precision> {
        if dist > self.max_dist {
            return None;
        }
        Some(self.max_dist - dist)
    }

    fn at_center(&self) -> Precision {
        self.max_dist
    }
}

/// Kernel for Gaussian smoothing.
#[derive(Copy, Clone, Debug)]
pub struct Gaussian {
    double_variance: Precision,
    cut_off_weight: Precision,
    at_center: Precision,
}

impl Gaussian {
    pub fn new(stdev: Precision, width: Precision) -> Self {
        let variance = stdev * stdev;
        let cut_off_weight = gaussian_dist(variance, stdev * width);
        Self {
            double_variance: variance * 2.0,
            cut_off_weight,
            at_center: gaussian_dist2(variance, 0.0),
        }
    }
}

impl Kernel for Gaussian {
    fn at_center(&self) -> Precision {
        self.at_center
    }

    fn weigh_dist(&self, dist: Precision) -> Option<Precision> {
        self.weigh_dist2(dist * dist)
    }

    fn weigh_dist2(&self, dist2: Precision) -> Option<Precision> {
        let w = (-dist2 / self.double_variance).exp();
        if w > self.cut_off_weight {
            None
        } else {
            Some(w)
        }
    }
}

fn gaussian_dist2(variance: Precision, dist2: Precision) -> Precision {
    // skip the constant term (peak height); it gets normalised out in the weighted mean
    // use dist2 because it's faster to calculate, and it saves us having to square it again here
    // mean = 0, val will always be positive
    (-dist2 / (2.0 * variance)).exp()
}

fn gaussian_dist(variance: Precision, dist: Precision) -> Precision {
    gaussian_dist2(variance, dist * dist)
}

struct WeightCache<'a, K: Kernel, const D: usize> {
    line: &'a [Point<Precision, D>],
    kernel: K,
    // would be lower memory if this were LRU
    cache: HashMap<(usize, usize), Option<Precision>>,
}

impl<'a, K: Kernel, const D: usize> WeightCache<'a, K, D> {
    pub fn new(line: &'a [Point<Precision, D>], kernel: K) -> Self {
        Self {
            line,
            kernel,
            cache: HashMap::default(),
        }
    }

    /// None if the weight of the edge is too low.
    /// Does not enforce key order or check that indices are in range.
    pub fn get_weight_unchecked(&mut self, idx1: usize, idx2: usize) -> Option<Precision> {
        let key = (idx1, idx2);
        let w;
        if let Entry::Vacant(e) = self.cache.entry(key) {
            w = self
                .kernel
                .weigh_dist2(distance_squared(&self.line[idx1], &self.line[idx2]));
            e.insert(w);
        } else {
            w = *self.cache.get(&key).unwrap();
        }
        w
    }

    /// None if the indices are out of range,
    /// or if the weight of the edge is too low.
    /// Enforces key order.
    pub fn get_weight(&mut self, idx1: usize, idx2: usize) -> Option<Precision> {
        let (lesser, greater) = match idx1.cmp(&idx2) {
            Ordering::Less => (idx1, idx2),
            Ordering::Equal => return Some(self.at_center()),
            Ordering::Greater => (idx2, idx1),
        };
        if lesser >= self.line.len() {
            return None;
        }
        self.get_weight_unchecked(lesser, greater)
    }

    pub fn at_center(&self) -> Precision {
        self.kernel.at_center()
    }
}

fn reflect_point<const D: usize>(
    reflect: &Point<Precision, D>,
    reflect_around: &Point<Precision, D>,
) -> Point<Precision, D> {
    2.0 * reflect_around - reflect.coords
}

/// Smooth line by applying an arbitrary kernel.
pub fn smooth_convolve<K: Kernel, const D: usize>(
    line: &[Point<Precision, D>],
    kernel: K,
) -> Vec<Point<Precision, D>> {
    let mut weight_cache = WeightCache::new(line, kernel);

    let first_point = line.first().unwrap();
    let last_point = line.last().unwrap();
    let last_idx = line.len() - 1;

    // We want the end points to stay where they are.
    // If we just cut off the kernel near the ends of the line,
    // the points would be smoothed away from the end points.
    // So we create some fake points just before and after the line,
    // effectively reflections of the first few and last few points,
    // to balance out the smoothing.

    let reflected_l: Vec<_> = (1..)
        .map(|idx| {
            weight_cache
                .get_weight(0, idx)
                .map(|w| (reflect_point(&line[idx], first_point), w))
        })
        .take_while(|o| o.is_some())
        .map(|o| o.unwrap())
        .collect();

    let reflected_r: Vec<_> = (1..)
        .map(|idx| {
            weight_cache
                .get_weight(last_idx, last_idx - idx)
                .map(|w| (reflect_point(&line[last_idx - idx], last_point), w))
        })
        .take_while(|o| o.is_some())
        .map(|o| o.unwrap())
        .collect();

    let mut smoothed = vec![];

    for (current_idx, current_point) in line.iter().enumerate() {
        let mut these_points = vec![(*current_point, weight_cache.at_center())];
        let mut to_reflect: usize = 0;

        // Go forward from the current point, possibly off the end of the line
        for idx_diff in 1.. {
            let next_idx = current_idx + idx_diff;

            if next_idx >= line.len() {
                // points will need to be taken from the pre-generated reflected points
                to_reflect = next_idx - line.len();
                break;
            }

            // can use *_unchecked because we know the index order and that they're in range
            if let Some(weight) = weight_cache.get_weight_unchecked(current_idx, next_idx) {
                these_points.push((line[next_idx], weight));
            } else {
                break;
            }
        }

        these_points.extend(reflected_r.iter().take(to_reflect));

        to_reflect = 0;

        // Go backward from the current point, possibly off the start of the line
        for idx_diff in 1.. {
            if current_idx < idx_diff {
                to_reflect = idx_diff - current_idx;
                break;
            }
            let next_idx = current_idx - idx_diff;

            if let Some(weight) = weight_cache.get_weight_unchecked(current_idx, next_idx) {
                these_points.push((line[next_idx], weight));
            } else {
                break;
            }
        }

        these_points.extend(reflected_l.iter().take(to_reflect));

        smoothed.push(weighted_mean(&these_points[..]));
    }

    smoothed
}
