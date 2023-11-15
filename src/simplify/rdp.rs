//! Simplify a linestring using the [Ramer-Douglas-Peucker](https://en.wikipedia.org/wiki/Ramer%E2%80%93Douglas%E2%80%93Peucker_algorithm) algorithm.
use crate::Precision;
use nalgebra::{distance_squared, Point};

fn proj_dist2<const D: usize>(
    start: &Point<Precision, D>,
    end: &Point<Precision, D>,
    p: &Point<Precision, D>,
    length_sq: Precision,
) -> Precision {
    let lensq_along = (p - start).dot(&(end - start));
    let proj;
    if lensq_along <= 0.0 {
        proj = start;
    } else if lensq_along >= length_sq {
        proj = end;
    } else {
        return distance_squared(p, &(start + lensq_along / length_sq * (end - start)));
    }
    distance_squared(p, proj)
}

fn rdp_keep_inner<const D: usize>(
    line: &[Point<Precision, D>],
    epsilon_sq: Precision,
    offset: usize,
) -> Vec<usize> {
    if line.len() <= 2 {
        return vec![];
    }

    let first = line.first().unwrap();
    let last = line.last().unwrap();
    let length_sq = distance_squared(first, last);

    let mut greatest_dist2 = (0, Precision::NEG_INFINITY);
    for (idx, point) in line.iter().enumerate().skip(1).take(line.len() - 2) {
        let d2 = proj_dist2(first, last, point, length_sq);
        if d2 > greatest_dist2.1 {
            greatest_dist2 = (idx + offset, d2)
        }
    }

    let mut to_keep = vec![];
    if greatest_dist2.1 > epsilon_sq {
        if greatest_dist2.0 > 1 {
            to_keep.append(&mut rdp_keep_inner(
                &line[0..=greatest_dist2.0],
                epsilon_sq,
                offset,
            ));
        }
        to_keep.push(greatest_dist2.0);
        if greatest_dist2.0 < line.len() - 2 {
            to_keep.append(&mut rdp_keep_inner(
                &line[greatest_dist2.0..line.len()],
                epsilon_sq,
                offset + greatest_dist2.0,
            ));
        }
    }

    to_keep
}

/// Return the indices of the points in the line which would be kept if simplified using RDP.
pub fn rdp_keep<const D: usize>(line: &[Point<Precision, D>], epsilon: Precision) -> Vec<usize> {
    let epsilon_sq = epsilon * epsilon;
    let mut out = Vec::with_capacity(line.len());
    out.push(0);
    out.append(&mut rdp_keep_inner(line, epsilon_sq, 0));
    out.push(line.len() - 1);
    out
}

/// Decimate the linestring using RDP.
pub fn rdp_reduce<const D: usize>(
    line: &[Point<Precision, D>],
    epsilon: Precision,
) -> Vec<Point<Precision, D>> {
    let kept = rdp_keep(line, epsilon);
    kept.into_iter().map(|idx| line[idx]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::test_utils::make_line;

    fn assert_reduce(orig: Vec<[f64; 2]>, expected: Vec<[f64; 2]>, epsilon: f64) {
        let orig_line = make_line(orig);
        let exp_line = make_line(expected);

        let out = rdp_reduce(orig_line.as_slice(), epsilon);
        assert_eq!(out, exp_line);
    }

    #[test]
    fn reduce() {
        assert_reduce(
            vec![[0.0, 0.0], [1.0, 0.1], [2.0, 0.0]],
            vec![[0.0, 0.0], [2.0, 0.0]],
            0.2,
        );
    }

    #[test]
    fn reduce_multi() {
        assert_reduce(
            vec![[0.0, 0.0], [0.5, 0.6], [1.0, 1.0], [1.6, 0.5], [2.0, 0.0]],
            vec![[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]],
            0.2,
        )
    }
}
