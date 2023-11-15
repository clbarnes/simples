use lru::LruCache;
use nalgebra::{distance, Point};

use crate::Precision;

/// Find length of a segment between any two points of a linestring.
///
/// Optionally caches pairwise distances.
struct DistanceFinder<'a, const D: usize> {
    points: &'a [Point<Precision, D>],
    closed: bool,
    cache: Option<LruCache<(usize, usize), Precision>>,
    total_length: Option<Precision>,
}

impl<'a, const D: usize> DistanceFinder<'a, D> {
    pub fn new(points: &'a [Point<Precision, D>], closed: bool, cache_size: Option<usize>) -> Self {
        Self {
            points,
            closed,
            cache: cache_size.map(|s| LruCache::new(NonZeroUsize::new(s).unwrap())),
            total_length: None,
        }
    }

    pub fn total_length(&mut self) -> Precision {
        if let Some(t) = self.total_length {
            return t;
        }

        if self.points.len() < 2 {
            return 0.0;
        }
        let mut out = (0..self.points.len() - 1)
            .map(|a| self.length_pair(a, a + 1))
            .sum::<Precision>();
        if self.closed {
            out += self.length_pair(self.points.len() - 1, 0);
        }
        self.total_length = Some(out);

        out
    }

    /// Ascending order.
    ///
    /// For closed linestrings, the closing pair is (last_idx, 0)
    fn length_pair(&mut self, start: usize, end: usize) -> Precision {
        if let Some(c) = self.cache.as_mut() {
            let key = (start, end);
            let r = c.get_or_insert(key, || distance(&self.points[start], &self.points[end]));
            *r
        } else {
            distance(&self.points[start], &self.points[end])
        }
    }

    /// Length in ascending order, i.e. `start` must be lower than `end` on a non-closed linestring.
    /// [None] if indices are not present or `end < start` and the linestring is not closed.
    pub fn length(&mut self, start: usize, end: usize) -> Option<Precision> {
        use std::cmp::Ordering::*;
        if start >= self.points.len() || end >= self.points.len() {
            return None;
        }

        match start.cmp(&end) {
            Less => Some(
                (start..end)
                    .map(|i| self.length_pair(i, i + 1))
                    .sum::<Precision>(),
            ),
            Equal => Some(0.0),
            Greater => {
                if !self.closed {
                    return None;
                }
                let out = (start..(self.points.len() - 1))
                    .chain(0..end)
                    .map(|i| self.length_pair(i, i + 1))
                    .sum::<Precision>();
                Some(out + self.length_pair(self.points.len() - 1, 0))
            }
        }
    }
}
