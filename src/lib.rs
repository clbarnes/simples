//! Simplify and smooth linestrings in N dimensions.
pub use nalgebra;
pub use nalgebra::Point;
pub use num_traits::Float;

pub mod simplify;
pub mod smooth;

pub type Precision = f64;

pub type Coord<T, const D: usize> = [T; D];

pub trait Location<T: Float, const D: usize> {
    fn location(&self) -> Coord<T, D>;

    /// The squared euclidean distance to another location.
    fn distance2_to<L: Location<T, D>>(&self, other: &L) -> T {
        sum(self
            .location()
            .iter()
            .zip(other.location().iter())
            .map(|(a, b)| (*a - *b).powi(2)))
    }

    /// The euclidean distance to another location.
    fn distance_to<L: Location<T, D>>(&self, other: &L) -> T {
        self.distance2_to(other).sqrt()
    }

    /// Where you would end up if you travelled `distance` towards `other`,
    /// and the overshoot: how far past the point you have travelled
    /// (negative if the point was not reached).
    fn project_towards<L: Location<T, D>>(&self, other: &L, distance: T) -> (Coord<T, D>, T) {
        let self_loc = self.location();
        let distance_to = self.distance_to(other);
        if (distance_to * distance).is_zero() {
            return (self_loc, T::zero());
        }
        let out = zip_op(&self_loc, &other.location(), |a, b| {
            let diff = *b - *a;
            *a + (diff / distance_to) * distance
        });
        (out, distance - distance_to)
    }
}

pub trait LocationMut<T: Float, const D: usize>: Location<T, D> {
    fn set_location(&mut self, new_loc: [T; D]);
}

impl<T: Float, const D: usize> Location<T, D> for Coord<T, D> {
    fn location(&self) -> Coord<T, D> {
        *self
    }
}

impl<T: Float, const D: usize> LocationMut<T, D> for Coord<T, D> {
    fn set_location(&mut self, new_loc: [T; D]) {
        *self = new_loc;
    }
}

impl<T: Float, const D: usize> Location<T, D> for &Coord<T, D> {
    fn location(&self) -> Coord<T, D> {
        **self
    }
}

impl<T: Float, const D: usize> Location<T, D> for &mut Coord<T, D> {
    fn location(&self) -> Coord<T, D> {
        **self
    }
}

impl<T: Float, const D: usize> LocationMut<T, D> for &mut Coord<T, D> {
    fn set_location(&mut self, new_loc: [T; D]) {
        **self = new_loc;
    }
}

fn sum<T: Float, I: IntoIterator<Item = T>>(it: I) -> T {
    it.into_iter().fold(T::zero(), |tot, this| tot + this)
}

fn zip_op<T: Float, const D: usize, F: Fn(&T, &T) -> T>(
    lhs: &Coord<T, D>,
    rhs: &Coord<T, D>,
    f: F,
) -> Coord<T, D> {
    let mut out = [T::zero(); D];
    lhs.iter()
        .zip(rhs.iter())
        .enumerate()
        .for_each(|(idx, (l, r))| {
            out[idx] = f(l, r);
        });
    out
}

fn sqdist<T: Float, I: IntoIterator<Item = T>>(v: I) -> T {
    sum(v.into_iter().map(|x| x.powi(2)))
}

fn dist<T: Float, I: IntoIterator<Item = T>>(v: I) -> T {
    sqdist(v).sqrt()
}

/// Find the total length of a linestring.
pub fn total_length<T: Float, const D: usize>(line: &[Coord<T, D>]) -> T {
    if line.len() < 2 {
        return T::zero();
    }
    sum(line
        .windows(2)
        .map(|points| points[0].distance_to(&points[1])))
}

#[cfg(test)]
mod test_utils {
    use nalgebra::Point2;
    pub type Pt = Point2<f64>;

    pub fn make_line(arrs: Vec<[f64; 2]>) -> Vec<Pt> {
        arrs.into_iter().map(|p| p.into()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn length1() {
        let line: Vec<_> = (0..5).map(|n| [n as f64]).collect();
        assert_eq!(total_length(&line), 4.0);
    }

    #[test]
    fn length2() {
        let line: Vec<_> = (0..2).map(|n| [n as f64; 2]).collect();
        assert_eq!(total_length(&line), 2.0f64.sqrt());
    }

    #[test]
    fn length3() {
        let line: Vec<_> = (0..2).map(|n| [n as f64; 3]).collect();
        assert_eq!(total_length(&line), 3.0f64.sqrt());
    }
}
