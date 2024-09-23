use lib::field::{print_elliptic_curve_points, Point};

fn main() {
    let g = Point::Point { x: 1, y: 2 };

    print_elliptic_curve_points(g);
}
