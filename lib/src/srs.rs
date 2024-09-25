use crate::{
    field::{point_add, scalar_mult, Point, F101},
    field_extension::F1012,
    polynomial::Polynomial,
};

pub struct SRS {
    pub f101_points: Vec<Point>,
    pub f101_2_elements: Vec<F1012>,
}

pub fn get_srs() -> SRS {
    let g1 = Point::Point {
        x: F101::new(1),
        y: F101::new(2),
    };
    let mut srs = Vec::new();

    // compute 1*G1, 2*G1, 2^2*G1, 2^3*G1,..., 2^6*G1
    // tau^0*G1, tau^1*G1...
    for i in 0..=6 {
        let k = 2u32.pow(i);
        let point = scalar_mult(k, g1);
        // println!("2^{} * G1 = {:?}", i, point);
        srs.push(point);
    }
    // two F101_2 element here for pairing, G2, tau*G2 (2*G2)
    let g2 = F1012::new(36, 31);
    let g2_2 = F1012::new(90, 82);

    SRS {
        f101_points: srs,
        f101_2_elements: [g2, g2_2].to_vec(),
    }
}

pub fn commit_poly(poly_for_commit: Polynomial) -> Point {
    let srs = get_srs();
    let mut committed_point = Point::Infinity;
    for j in 0..poly_for_commit.coeffs.len() {
        let coeff = poly_for_commit.coeffs[j].0;
        let p = scalar_mult(coeff, srs.f101_points[j]);
        committed_point = point_add(committed_point, p);
    }
    committed_point
}
