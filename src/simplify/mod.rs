//! Simplify linestrings.
//!
//! Linestrings are generally simplified by changing the number of points.
pub mod rdp;
pub mod vw;

use crate::{total_length, Precision};
use nalgebra::Point;
use std::cmp::{Ordering, PartialOrd};

/// Traverse the edges of a linestring, dropping a new node every `sample_distance`, and deleting all original nodes.
/// `offset` allows you to start partway down the first edge: use 0.0 if you want to include the first node.
pub fn sample_every<const D: usize>(
    line: &[Point<Precision, D>],
    sample_distance: Precision,
    offset: Precision,
) -> Vec<Point<Precision, D>> {
    let mut iter = line.iter();
    if line.len() <= 1 {
        return iter.cloned().collect();
    }
    let mut prev = *iter.next().unwrap();
    let mut out = Vec::default();
    let mut remaining_dist: f64;
    if offset == 0.0 {
        out.push(prev);
        remaining_dist = offset;
    } else {
        remaining_dist = sample_distance
    }
    let mut next = *iter.next().unwrap();

    loop {
        let vec = next - prev;
        let edge_length = vec.magnitude();

        remaining_dist = match remaining_dist.partial_cmp(&edge_length).unwrap() {
            Ordering::Less => {
                prev += (vec / edge_length) * remaining_dist;
                out.push(prev);
                sample_distance
            }
            Ordering::Equal => {
                prev = next;
                out.push(prev);
                let next_opt = iter.next();
                if next_opt.is_none() {
                    break;
                }
                next = *next_opt.unwrap();
                sample_distance
            }
            Ordering::Greater => {
                prev = next;
                let next_opt = iter.next();
                if next_opt.is_none() {
                    out.push(prev);
                    break;
                }
                next = *next_opt.unwrap();
                remaining_dist - edge_length
            }
        };
    }

    out
}

/// Resample a linestring to ensure that it has `n_points` points,
/// by dividing the total length evenly.
pub fn resample<const D: usize>(
    line: &[Point<Precision, D>],
    n_points: usize,
) -> Vec<Point<Precision, D>> {
    let len = total_length(line);
    if len == 0.0 {
        panic!("Not enough points");
    }
    let dist = len / ((n_points - 1) as f64);
    sample_every(line, dist, 0.0)
}
