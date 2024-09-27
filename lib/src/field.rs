use std::{collections::HashMap, ops::Sub};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct F101(pub u32);

impl F101 {
    pub const P: u32 = 101;

    pub const ZERO: F101 = F101(0);
    pub const ONE: F101 = F101(1);
    pub const NEG_ONE: F101 = F101(100);

    pub fn new(n: u32) -> Self {
        F101(n % Self::P)
    }

    pub fn add(self, other: F101) -> F101 {
        F101::new((self.0 + other.0) % Self::P)
    }

    pub fn sub(self, other: F101) -> F101 {
        F101::new((self.0 + Self::P - other.0) % Self::P)
    }

    pub fn mul(self, other: F101) -> F101 {
        F101::new((self.0 * other.0) % Self::P)
    }

    pub fn inv(self) -> Option<F101> {
        if self.0 == 0 {
            None
        } else {
            Some(F101::new(mod_pow(self.0, Self::P - 2, Self::P)))
        }
    }

    pub fn value(self) -> u32 {
        self.0
    }

    pub fn neg(self) -> F101 {
        F101::new(F101::P - self.0)
    }
}

pub fn mod_pow(mut base: u32, mut exp: u32, modulus: u32) -> u32 {
    let mut result = 1u32;
    base %= modulus;
    while exp > 0 {
        if exp % 2 == 1 {
            result = (result * base) % modulus;
        }
        exp >>= 1;
        base = (base * base) % modulus;
    }
    result
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Point {
    Infinity,
    Point { x: F101, y: F101 },
}

impl Point {
    pub fn point_neg(self) -> Point {
        match self {
            Point::Infinity => Point::Infinity,
            Point::Point { x, y } => Point::Point { x, y: y.neg() },
        }
    }

    pub fn add(self, other: Point) -> Point {
        point_add(self, other)
    }

    pub fn mul(self, k: F17) -> Point {
        scalar_mult(k.0, self)
    }
}

pub fn point_add(p1: Point, p2: Point) -> Point {
    match (p1, p2) {
        (Point::Infinity, _) => p2,
        (_, Point::Infinity) => p1,
        (Point::Point { x: x1, y: y1 }, Point::Point { x: x2, y: y2 }) => {
            if x1 == x2 && (y1.add(y2)).value() == 0 {
                Point::Infinity
            } else {
                let m = if x1 == x2 && y1 == y2 {
                    let numerator = x1.mul(x1).mul(F101::new(3));
                    let denominator = y1.mul(F101::new(2));
                    match denominator.inv() {
                        Some(inv) => numerator.mul(inv),
                        None => return Point::Infinity,
                    }
                } else {
                    let numerator = y2.sub(y1);
                    let denominator = x2.sub(x1);
                    match denominator.inv() {
                        Some(inv) => numerator.mul(inv),
                        None => return Point::Infinity,
                    }
                };
                let x3 = m.mul(m).sub(x1).sub(x2);
                let y3 = m.mul(x1.sub(x3)).sub(y1);
                Point::Point { x: x3, y: y3 }
            }
        }
    }
}

pub fn scalar_mult(mut k: u32, p: Point) -> Point {
    let mut result = Point::Infinity;
    let mut addend = p;

    while k > 0 {
        if k % 2 == 1 {
            result = point_add(result, addend);
        }
        addend = point_add(addend, addend);
        k >>= 1;
    }
    result
}

pub fn print_elliptic_curve_points(g: Point) {
    let mut points = HashMap::new();
    let mut order = 1u32;
    let mut p = g;

    loop {
        points.insert(order, p);
        order += 1;
        p = point_add(p, g);
        if p == Point::Infinity {
            points.insert(order, Point::Infinity);
            break;
        }
    }

    println!(
        "There are in total {} points on the elliptic curve (including infinite point).",
        points.len()
    );

    for k in 1..=order {
        if let Some(point) = points.get(&k) {
            match point {
                Point::Infinity => println!("{}: Infinity", k),
                Point::Point { x, y } => println!("{}: ({}, {})", k, x.value(), y.value()),
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct F17(pub u32);

impl F17 {
    pub const P: u32 = 17;

    pub const ZERO: F17 = F17(0);
    pub const ONE: F17 = F17(1);
    pub const NEG_ONE: F17 = F17(16);

    pub const H: [F17; 4] = [F17(1), F17(4), F17(16), F17(13)];
    pub const K1H: [F17; 4] = [F17(2), F17(8), F17(15), F17(9)];
    pub const K2H: [F17; 4] = [F17(3), F17(12), F17(14), F17(5)];

    pub const K1: F17 = F17(2);
    pub const K2: F17 = F17(3);

    pub fn new(n: u32) -> Self {
        F17(n % Self::P)
    }

    pub fn add(self, other: F17) -> F17 {
        F17::new((self.0 + other.0) % Self::P)
    }

    pub fn sub(self, other: F17) -> F17 {
        F17::new((self.0 + Self::P - other.0) % Self::P)
    }

    pub fn mul(self, other: F17) -> F17 {
        F17::new((self.0 * other.0) % Self::P)
    }

    pub fn inv(self) -> Option<F17> {
        if self.0 == 0 {
            None
        } else {
            Some(F17::new(mod_pow(self.0, Self::P - 2, Self::P)))
        }
    }

    pub fn value(self) -> u32 {
        self.0
    }

    pub fn pow(self, exp: u32) -> F17 {
        let mut result = F17::ONE;
        for _ in 0..exp {
            result = result.mul(self);
        }
        result
    }

    pub fn neg(self) -> F17 {
        F17::ZERO.sub(self)
    }
}

// fn vander_matrix(points: &[(F17, F17)]) -> [[F17; 4]; 4] {
//     let mut matrix = [[F17::new(0); 4]; 4];
//     for i in 0..4 {
//         let x = points[i].0;
//         matrix[i][0] = x.mul(x).mul(x); // x^3
//         matrix[i][1] = x.mul(x); // x^2
//         matrix[i][2] = x; // x^1
//         matrix[i][3] = F17::new(1); // x^0
//     }
//     matrix
// }

fn vander_matrix(points: &[(F17, F17)]) -> [[F17; 4]; 4] {
    let mut matrix = [[F17::new(0); 4]; 4];
    for i in 0..4 {
        let x = points[i].0;
        matrix[i][0] = F17::new(1); // x^0
        matrix[i][1] = x; // x^1
        matrix[i][2] = x.mul(x); // x^2
        matrix[i][3] = x.mul(x).mul(x); // x^3
    }
    matrix
}

fn matrix_inverse(matrix: [[F17; 4]; 4]) -> Option<[[F17; 4]; 4]> {
    let mut augmented_matrix = [[F17::new(0); 8]; 4];

    // [matrix | identity] augmented_matrix
    for i in 0..4 {
        for j in 0..4 {
            augmented_matrix[i][j] = matrix[i][j]; // left side: original matrix
        }
        augmented_matrix[i][i + 4] = F17::new(1); // right side: unit matrix
    }

    // gauss elimination
    for i in 0..4 {
        // find pivot (non-zero)
        if augmented_matrix[i][i] == F17::new(0) {
            // If the pivot is zero, swap with a non-zero row.
            let mut swap_row = None;
            for j in i + 1..4 {
                if augmented_matrix[j][i] != F17::new(0) {
                    swap_row = Some(j);
                    break;
                }
            }
            if let Some(row) = swap_row {
                augmented_matrix.swap(i, row);
            } else {
                return None;
            }
        }

        // Change the pivot element to 1
        let inv_pivot = augmented_matrix[i][i].inv()?;
        for j in 0..8 {
            augmented_matrix[i][j] = augmented_matrix[i][j].mul(inv_pivot);
        }

        // Eliminate the other rows so that the remaining elements of the column are 0
        for j in 0..4 {
            if i != j {
                let factor = augmented_matrix[j][i];
                for k in 0..8 {
                    augmented_matrix[j][k] =
                        augmented_matrix[j][k].sub(augmented_matrix[i][k].mul(factor));
                }
            }
        }
    }

    // the inverse matrix on the right side
    let mut inverse_matrix = [[F17::new(0); 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            inverse_matrix[i][j] = augmented_matrix[i][j + 4];
        }
    }

    Some(inverse_matrix)
}

fn matrix_mul_vec(matrix: [[F17; 4]; 4], vec: [F17; 4]) -> [F17; 4] {
    let mut result = [F17::new(0); 4];
    for i in 0..4 {
        for j in 0..4 {
            result[i] = result[i].add(matrix[i][j].mul(vec[j]));
        }
    }
    result
}

// vandermonde_inverse * vec[points.y]
pub fn solve_coefficients(points: &[(F17, F17); 4]) -> Option<[F17; 4]> {
    let vandermonde = vander_matrix(points);
    // println!("vandermonde: {:?}", vandermonde);
    let inverse = matrix_inverse(vandermonde)?;
    // println!("vandermonde_inverse: {:?}", inverse);
    // let y_values: [F17; 4] = [F17::new(3), F17::new(4), F17::new(5), F17::new(9)];
    let y_values: [F17; 4] = points
        .into_iter()
        .map(|point| point.1)
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    let coefficients = matrix_mul_vec(inverse, y_values);
    Some(coefficients)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn generate_coset(h: &[F17], k: F17) -> Vec<F17> {
        h.iter().map(|&x| (x.mul(k))).collect()
    }

    #[test]
    fn coset_test() {
        println!("{:?}", generate_coset(&F17::H, F17::new(2)));
        println!("{:?}", generate_coset(&F17::H, F17::new(3)));

        assert_eq!(F17::K1H.to_vec(), generate_coset(&F17::H, F17::new(2)));
        assert_eq!(F17::K2H.to_vec(), generate_coset(&F17::H, F17::new(3)));
    }

    #[test]
    fn solve_coefficients_test() {
        let points: [(F17, F17); 4] = [
            (F17::new(1), F17::new(3)),
            (F17::new(4), F17::new(4)),
            (F17::new(16), F17::new(5)),
            (F17::new(13), F17::new(9)),
        ];
        let result = solve_coefficients(&points).unwrap();
        println!("result poly in coeff form: {:?}", result);
    }

    #[test]
    fn larange_base_test() {
        let points: [(F17, F17); 4] = [
            (F17::new(1), F17::new(1)),
            (F17::new(4), F17::new(0)),
            (F17::new(16), F17::new(0)),
            (F17::new(13), F17::new(0)),
        ];
        let result = solve_coefficients(&points).unwrap();
        println!("result poly in coeff form: {:?}", result);
    }
}
