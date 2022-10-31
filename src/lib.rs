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
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
