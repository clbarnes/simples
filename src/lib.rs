//! Simplify and smooth linestrings in N dimensions.
pub use nalgebra;
use nalgebra::{distance, Point};

pub mod simplify;
pub mod smooth;

pub type Precision = f64;

/// Find the total length of a linestring.
pub fn total_length<const D: usize>(line: &[Point<Precision, D>]) -> Precision {
    if line.len() < 2 {
        return 0.0;
    }
    line.windows(2)
        .map(|points| distance(&points[0], &points[1]))
        .sum::<Precision>()
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Point1, Point2, Point3};

    #[test]
    fn length1() {
        let line: Vec<_> = (0..5).map(|n| Point1::new(n as f64)).collect();
        assert_eq!(total_length(&line), 4.0);
    }

    #[test]
    fn length2() {
        let line: Vec<_> = (0..2).map(|n| Point2::new(n as f64, n as f64)).collect();
        assert_eq!(total_length(&line), (2.0 as f64).sqrt());
    }

    #[test]
    fn length3() {
        let line: Vec<_> = (0..2)
            .map(|n| Point3::new(n as f64, n as f64, n as f64))
            .collect();
        assert_eq!(total_length(&line), (3.0 as f64).sqrt());
    }
}
