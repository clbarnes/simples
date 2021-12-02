use nalgebra::{distance, distance_squared, Point};
use std::cmp::{Ord, Ordering, PartialOrd};
use std::collections::{BinaryHeap, HashSet};

pub type Precision = f64;

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

    let mut greatest_dist2 = (0, -Precision::INFINITY);
    for (idx, point) in line.iter().enumerate().skip(1).take(line.len() - 2) {
        let d2 = proj_dist2(first, last, point, length_sq);
        if d2 > greatest_dist2.1 {
            greatest_dist2 = (idx, d2)
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

pub fn rdp_keep<const D: usize>(line: &[Point<Precision, D>], epsilon: Precision) -> Vec<usize> {
    let epsilon_sq = epsilon * epsilon;
    let mut out = vec![0];
    out.append(&mut rdp_keep_inner(line, epsilon_sq, 0));
    out.push(line.len() - 1);
    out
}

pub fn rdp_reduce<const D: usize>(
    line: &[Point<Precision, D>],
    epsilon: Precision,
) -> Vec<Point<Precision, D>> {
    let kept = rdp_keep(line, epsilon);
    kept.into_iter().map(|idx| line[idx]).collect()
}

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

#[derive(Clone)]
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

    fn is_valid(&self, skipped: &HashSet<usize>) -> bool {
        let is_invalid = skipped.contains(&self.indices.0)
            || skipped.contains(&self.indices.1)
            || skipped.contains(&self.indices.2);
        !is_invalid
    }

    fn get_replacement(&self, points: &[Point<Precision, D>], skipped: &HashSet<usize>) -> Self {
        let (new_left, new_right) = neighbours_not_in(self.indices.0, self.indices.2, skipped);
        Self::from_indices(points, (new_left, self.indices.1, new_right))
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
        other.area.partial_cmp(&self.area)
    }
}

impl<const D: usize> Ord for Triangle<D> {
    fn cmp(&self, other: &Self) -> Ordering {
        // operands reversed so it's a min-heap
        other.partial_cmp(other).unwrap_or(Ordering::Equal)
    }
}

fn neighbours_not_in(
    starting_left: usize,
    starting_right: usize,
    skipped: &HashSet<usize>,
) -> (usize, usize) {
    let mut left = starting_left;
    while skipped.contains(&left) {
        left -= 1;
    }
    let mut right = starting_right;
    while skipped.contains(&right) {
        right += 1;
    }
    (left, right)
}

fn vw_drop<const D: usize>(line: &[Point<Precision, D>], n_points: usize) -> HashSet<usize> {
    if line.len() <= 2.max(n_points) || line.len() >= n_points {
        return HashSet::default();
    }
    let mut queue = BinaryHeap::default();
    let mut drop = HashSet::default();
    for idx in 0..(line.len() - 2) {
        queue.push(Triangle::from_indices(line, (idx, idx + 1, idx + 2)))
    }
    while line.len() - drop.len() > n_points {
        if let Some(tri) = queue.pop() {
            if tri.is_valid(&drop) {
                drop.insert(tri.center_index());
            } else {
                queue.push(tri.get_replacement(line, &drop));
            }
        } else {
            break;
        }
    }
    drop
}

pub fn vw_keep<const D: usize>(line: &[Point<Precision, D>], n_points: usize) -> Vec<usize> {
    let drop = vw_drop(line, n_points);
    (0..line.len()).filter(|idx| !drop.contains(idx)).collect()
}

pub fn vw_reduce<const D: usize>(
    line: &[Point<Precision, D>],
    n_points: usize,
) -> Vec<Point<Precision, D>> {
    let drop = vw_drop(line, n_points);
    line.iter()
        .enumerate()
        .filter_map(|(idx, p)| if !drop.contains(&idx) { None } else { Some(*p) })
        .collect()
}

pub fn sample_every<const D: usize>(
    line: &[Point<Precision, D>],
    sample_distance: Precision,
) -> Vec<Point<Precision, D>> {
    let mut iter = line.iter();
    if line.len() <= 1 {
        return iter.cloned().collect();
    }
    let mut prev = *iter.next().unwrap();
    let mut out = vec![prev];
    let mut next = *iter.next().unwrap();

    // how far until we drop another point
    let mut remaining_dist = sample_distance;

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

pub fn total_length<const D: usize>(line: &[Point<Precision, D>]) -> Precision {
    if line.len() < 2 {
        return 0.0;
    }
    line.windows(2)
        .map(|points| distance(&points[0], &points[1]))
        .sum::<Precision>()
}

pub fn resample<const D: usize>(
    line: &[Point<Precision, D>],
    n_points: usize,
) -> Vec<Point<Precision, D>> {
    let len = total_length(line);
    if len == 0.0 {
        panic!("Not enough points");
    }
    let dist = len / ((n_points - 1) as f64);
    sample_every(line, dist)
}
