//! Resample a linestring by placing evenly-spaced points along its length.
use crate::{total_length, Precision};
use nalgebra::Point;
use std::cmp::{Ordering, PartialOrd};

/// Create a new linestring by traversing the original, placing a node every `sample_distance`.
/// `offset` allows you to start partway down the first edge: use 0.0 if you want to include the first node.
/// Returns the resampled points and the geodesic distance from the last resampled point to the original last point.
///
/// A zero-point line remains zero-point; a single-point line keeps that single point.
/// A line shorter than `sample_distance + offset` will be reduced to a single point.
///
/// `sample_distance` must be positive and `offset` must be non-negative (panics if these are invalid).
pub fn sample_every<const D: usize>(
    line: &[Point<Precision, D>],
    sample_distance: Precision,
    offset: Precision,
) -> (Vec<Point<Precision, D>>, Precision) {
    if sample_distance <= 0.0 {
        panic!("`sample_distance` must be positive");
    }
    if offset < 0.0 {
        panic!("`offset` must be non-negative");
    }
    let mut iter = line.iter();
    if line.len() <= 1 {
        return (iter.cloned().collect(), 0.0);
    }
    let mut prev = *iter.next().unwrap();
    let mut out = Vec::default();
    let mut remaining_dist: f64;
    if offset == 0.0 {
        out.push(prev);
        remaining_dist = sample_distance;
    } else {
        remaining_dist = offset
    }
    let mut next = *iter.next().unwrap();

    loop {
        let vec = next - prev;
        let edge_length = vec.magnitude();

        match remaining_dist.partial_cmp(&edge_length).unwrap() {
            Ordering::Less => {
                prev += (vec / edge_length) * remaining_dist;
                out.push(prev);
                remaining_dist = sample_distance;
            }
            Ordering::Equal => {
                prev = next;
                out.push(prev);
                let Some(next_ref) = iter.next() else {
                    remaining_dist = 0.0;
                    break;
                };
                next = *next_ref;
                remaining_dist = sample_distance;
            }
            Ordering::Greater => {
                prev = next;
                remaining_dist -= edge_length;

                let Some(next_ref) = iter.next() else { break };
                next = *next_ref;
            }
        };
    }

    (out, sample_distance - remaining_dist)
}

/// Resample a linestring to ensure that it has `n_points` points,
/// by dividing the total length evenly.
///
/// Panics if line has zero length.
pub fn resample<const D: usize>(
    line: &[Point<Precision, D>],
    n_points: usize,
) -> Vec<Point<Precision, D>> {
    let len = total_length(line);
    if len == 0.0 {
        panic!("Not enough points");
    }
    let dist = len / ((n_points - 1) as f64);
    sample_every(line, dist, 0.0).0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn half_line() {
        let ls1: Vec<Point<f64, 1>> = vec![[0.0].into(), [1.0].into()];
        let resampled = resample(ls1.as_slice(), 3);
        println!("{:?}", resampled);
        assert_eq!(resampled.len(), 3);
        assert_eq!(resampled[0], ls1[0]);
        assert_eq!(resampled[1], [0.5].into());
        assert_eq!(resampled[2], ls1[1]);
    }

    #[test]
    fn resample_line() {
        let ls1: Vec<Point<f64, 1>> = vec![[0.0].into(), [3.0].into()];
        let (resampled, remainder) = sample_every(ls1.as_slice(), 1.0, 0.5);
        println!("{:?}", resampled);
        assert_eq!(resampled.len(), 3);
        assert_eq!(resampled[0], [0.5].into());
        assert_eq!(resampled[1], [1.5].into());
        assert_eq!(resampled[2], [2.5].into());
        assert_eq!(remainder, 0.5);
    }
}
