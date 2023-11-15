//! Simplify a linestring using the [Visvalingam-Whyatt](https://en.wikipedia.org/wiki/Visvalingam%E2%80%93Whyatt_algorithm) algorithm.
use crate::Precision;
use nalgebra::{distance, Point};
use std::cmp::{Ord, Ordering, PartialOrd};
use std::collections::{BinaryHeap, HashSet};

/// Area of a triangle based on the locations of its 3 corners.
fn tri_area<const D: usize>(
    p1: &Point<Precision, D>,
    p2: &Point<Precision, D>,
    p3: &Point<Precision, D>,
) -> Precision {
    let s1 = distance(p1, p2);
    let s2 = distance(p2, p3);
    let s3 = distance(p3, p1);

    let s = (s1 + s2 + s3) / 2.0;
    (s * (s - s1) * (s - s2) * (s - s3)).sqrt()
}

/// A triangle whose corners are indices into some reference table of points.
#[derive(Copy, Clone, Debug)]
struct Triangle<const D: usize> {
    pub indices: (usize, usize, usize),
    pub area: Precision,
}

impl<const D: usize> Triangle<D> {
    fn from_indices(all_points: &[Point<Precision, D>], indices: (usize, usize, usize)) -> Self {
        Self {
            indices,
            area: tri_area(
                &all_points[indices.0],
                &all_points[indices.1],
                &all_points[indices.2],
            ),
        }
    }

    /// Check whether the triangle is invalid, i.e. its corner indices appear in the given disallowed set.
    fn is_valid(&self, skipped: &HashSet<usize>) -> bool {
        let is_invalid = skipped.contains(&self.indices.0)
            || skipped.contains(&self.indices.1)
            || skipped.contains(&self.indices.2);
        !is_invalid
    }

    /// Find a triangle with the same midpoint and end points which do not appear in the `skipped` set.
    ///
    /// The left index goes leftwards, the right index goes rightwards.
    fn get_replacement(
        &self,
        points: &[Point<Precision, D>],
        skipped: &HashSet<usize>,
        wrapping: bool,
    ) -> Option<Self> {
        let (new_left, new_right) = neighbours_not_in(
            self.indices.0,
            self.indices.2,
            skipped,
            points.len(),
            wrapping,
        )?;
        Some(Self::from_indices(
            points,
            (new_left, self.indices.1, new_right),
        ))
    }

    fn center_index(&self) -> usize {
        self.indices.1
    }
}

impl<const D: usize> PartialEq for Triangle<D> {
    fn eq(&self, other: &Self) -> bool {
        self.area == other.area
    }
}

impl<const D: usize> Eq for Triangle<D> {}

impl<const D: usize> PartialOrd for Triangle<D> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // reversed operands for min-heap
        other.area.partial_cmp(&self.area)
    }
}

impl<const D: usize> Ord for Triangle<D> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap_or(Ordering::Equal)
    }
}

/// Scan left and right from `left` and `right` for non-`skipped` indices.
fn neighbours_not_in(
    mut left: usize,
    mut right: usize,
    skipped: &HashSet<usize>,
    len: usize,
    wrapping: bool,
) -> Option<(usize, usize)> {
    if skipped.len() >= len - 3 {
        return None;
    }
    while skipped.contains(&left) {
        if left > 0 {
            left -= 1;
        } else if wrapping {
            left = len - 1;
        } else {
            return None;
        }
    }
    while skipped.contains(&right) {
        if right == len - 1 {
            right += 1;
        } else if wrapping {
            right = 0;
        } else {
            return None;
        }
    }
    Some((left, right))
}

fn vw_drop<const D: usize>(
    line: &[Point<Precision, D>],
    n_points: usize,
    closed: bool,
) -> HashSet<usize> {
    if line.len() <= 2.min(n_points) {
        return HashSet::with_capacity(0);
    }
    let mut queue = BinaryHeap::default();
    let mut drop = HashSet::with_capacity(line.len() - n_points);
    for idx in 0..(line.len() - 2) {
        queue.push(Triangle::from_indices(line, (idx, idx + 1, idx + 2)))
    }
    if closed {
        let len = line.len();
        queue.push(Triangle::from_indices(line, (len - 2, len - 1, 0)));
        queue.push(Triangle::from_indices(line, (len - 1, 0, 1)));
    }
    while line.len() - drop.len() > n_points {
        if let Some(tri) = queue.pop() {
            if tri.is_valid(&drop) {
                drop.insert(tri.center_index());
            } else {
                if let Some(repl) = tri.get_replacement(line, &drop, closed) {
                    queue.push(repl);
                } else {
                    break;
                }
            }
        } else {
            break;
        }
    }
    drop
}

/// Return the indices of points on the linestring to be kept if decimated by VW.
///
/// `closed = true` where the linestring represents a polygon and there is an edge from the last point to the first.
pub fn vw_keep<const D: usize>(
    line: &[Point<Precision, D>],
    n_points: usize,
    closed: bool,
) -> Vec<usize> {
    let drop = vw_drop(line, n_points, closed);
    (0..line.len()).filter(|idx| !drop.contains(idx)).collect()
}

/// Decimate the linestring using VW.
///
/// `closed = true` where the linestring represents a polygon and there is an edge from the last point to the first.
pub fn vw_reduce<const D: usize>(
    line: &[Point<Precision, D>],
    n_points: usize,
    closed: bool,
) -> Vec<Point<Precision, D>> {
    let drop = vw_drop(line, n_points, closed);
    line.iter()
        .enumerate()
        .filter_map(|(idx, p)| if drop.contains(&idx) { None } else { Some(*p) })
        .collect()
}

#[cfg(test)]
mod tests {
    use nalgebra::Point2;

    type Pt = Point2<f64>;

    use super::*;

    fn make_line(arrs: Vec<[f64; 2]>) -> Vec<Pt> {
        arrs.into_iter().map(|p| p.into()).collect()
    }

    fn assert_reduce(
        orig: Vec<[f64; 2]>,
        expected: Vec<[f64; 2]>,
        n_points: Option<usize>,
        closed: bool,
    ) {
        let orig_line = make_line(orig);
        let exp_line = make_line(expected);

        let out = vw_reduce(
            orig_line.as_slice(),
            n_points.unwrap_or(exp_line.len()),
            closed,
        );
        assert_eq!(out, exp_line);
    }

    #[test]
    fn reduce() {
        assert_reduce(
            vec![[0.0, 0.0], [0.9, 0.0], [1.0, 1.0], [1.1, 0.0], [2.0, 0.0]],
            vec![[0.0, 0.0], [0.9, 0.0], [1.1, 0.0], [2.0, 0.0]],
            None,
            false,
        );
    }

    #[test]
    fn reduce_multi() {
        assert_reduce(
            vec![
                [0.0, 0.0],
                [0.9, 0.0],
                [1.0, 1.0],
                [1.1, 0.0],
                [1.9, 0.0],
                [2.0, 1.0],
                [2.1, 0.0],
                [3.0, 0.0],
            ],
            vec![
                [0.0, 0.0],
                [0.9, 0.0],
                [1.1, 0.0],
                [1.9, 0.0],
                [2.1, 0.0],
                [3.0, 0.0],
            ],
            None,
            false,
        )
    }

    #[test]
    fn reduce_wrap() {
        assert_reduce(
            vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [-0.1, 0.5]],
            vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
            None,
            true,
        );
    }

    #[test]
    fn reduce_nowrap() {
        assert_reduce(
            vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [-0.1, 0.5]],
            vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [-0.1, 0.5]],
            None,
            false,
        );
    }
}
