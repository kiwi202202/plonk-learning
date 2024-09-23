use std::collections::HashMap;

pub const P: u32 = 101;

// y^2 = x^3 + ax + b, a = 0, b = 3
pub const A: u32 = 0;
pub const B: u32 = 3;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Point {
    Infinity,
    Point { x: u32, y: u32 },
}

pub fn mod_inv(n: u32, p: u32) -> Option<u32> {
    if n == 0 {
        None
    } else {
        Some(mod_pow(n, p - 2, p))
    }
}

pub fn mod_pow(mut base: u32, mut exp: u32, modulus: u32) -> u32 {
    if modulus == 1 {
        return 0;
    }
    let mut result = 1u32;
    base = base % modulus;
    while exp > 0 {
        if exp % 2 == 1 {
            result = (result * base) % modulus;
        }
        exp = exp >> 1;
        base = (base * base) % modulus;
    }
    result
}

pub fn point_add(p1: Point, p2: Point) -> Point {
    match (p1, p2) {
        (Point::Infinity, _) => p2,
        (_, Point::Infinity) => p1,
        (Point::Point { x: x1, y: y1 }, Point::Point { x: x2, y: y2 }) => {
            if x1 == x2 && (y1 + y2) % P == 0 {
                Point::Infinity
            } else {
                let m = if x1 == x2 && y1 == y2 {
                    // double add
                    let numerator = (3 * x1 * x1 + A) % P;
                    let denominator = (2 * y1) % P;
                    match mod_inv(denominator, P) {
                        Some(inv) => (numerator * inv) % P,
                        None => return Point::Infinity,
                    }
                } else {
                    let numerator = (y2 + P - y1) % P;
                    let denominator = (x2 + P - x1) % P;
                    match mod_inv(denominator, P) {
                        Some(inv) => (numerator * inv) % P,
                        None => return Point::Infinity,
                    }
                };
                let x3 = (m * m + P - x1 + P - x2) % P;
                let y3 = (m * (x1 + P - x3) + P - y1) % P;
                Point::Point { x: x3, y: y3 }
            }
        }
    }
}

pub fn point_neg(p: Point) -> Point {
    match p {
        Point::Infinity => Point::Infinity,
        Point::Point { x, y } => Point::Point { x, y: (P - y) % P },
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
                Point::Point { x, y } => println!("{}: ({}, {})", k, x, y),
            }
        }
    }
}
