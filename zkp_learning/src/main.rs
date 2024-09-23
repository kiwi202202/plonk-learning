use lib::{
    field::{print_elliptic_curve_points, scalar_mult, Point, F101},
    field_extension::F1012,
};

fn main() {
    let g1 = Point::Point {
        x: F101::new(1),
        y: F101::new(2),
    };

    print_elliptic_curve_points(g1);

    let mut srs = Vec::new();

    // compute 1*G1, 2*G1, 2^2*G1, 2^3*G1,..., 2^6*G1
    // tau^0*G1, tau^1*G1...
    for i in 0..=6 {
        let k = 2u32.pow(i);
        let point = scalar_mult(k, g1);
        println!("2^{} * G1 = {:?}", i, point);
        srs.push(point);
    }

    // two F101_2 element here for pairing, G2, tau*G2 (2*G2)
    // srs.push(F1012::new(36, 31));
    // srs.push(F1012::new(90, 82));
}
