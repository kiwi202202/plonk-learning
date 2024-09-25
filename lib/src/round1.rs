use crate::{
    field::{Point, F17},
    polynomial::{get_Z_H, Polynomial},
    pythagorean_transcript::gen_t_polys_coeff,
    srs::commit_poly,
};

// generate random b1..b6 in F17
pub const B_RANDS: [F17; 6] = [F17(7), F17(4), F17(11), F17(12), F17(16), F17(2)];

pub fn compute_poly_coeff_round1(b_x: F17, b: F17, f: Polynomial) -> Polynomial {
    let zh_poly = get_Z_H();
    // compute (b1 * x + b2) * Z_H(x)
    let after_bx = zh_poly.mul_by_monomial(b_x, 1); // b_x * x * Z_H(x)
    let after_b_const = zh_poly.mul_by_monomial(b, 0); // b * Z_H(x)
    let final_zh = after_bx.add(&after_b_const); // b1 * x * Z_H(x) + b2 * Z_H(x)

    // a(x) = (b1x + b2) * Z_H(x) + f_a(x)
    final_zh.add(&f)
}

pub fn gen_round1_result() -> ([Polynomial; 3], [Point; 3]) {
    let fs = gen_t_polys_coeff();
    let mut round1_polys = Vec::new();
    let mut round1_committed_points = Vec::new();
    for i in 0..3 {
        let f = fs[i];
        let f_poly = Polynomial { coeffs: f.to_vec() };
        let poly_for_commit = compute_poly_coeff_round1(B_RANDS[2 * i], B_RANDS[2 * i + 1], f_poly);
        println!("{:?}", poly_for_commit);
        round1_polys.push(poly_for_commit.clone());
        // println!("Commited Point: {:?}", commit_poly(poly_for_commit));
        round1_committed_points.push(commit_poly(poly_for_commit));
    }
    (
        round1_polys.try_into().unwrap(),
        round1_committed_points.try_into().unwrap(),
    )
}

#[cfg(test)]
mod tests {

    use crate::{
        field::{point_add, scalar_mult, Point, F101},
        pythagorean_transcript::gen_t_polys_coeff,
        srs::{commit_poly, get_srs},
    };

    use super::*;

    // we commited the entire transcipt by the following 3 committed point, a-box, b-box, c-box
    #[test]
    fn compute_poly_coeff_test() {
        let fs = gen_t_polys_coeff();
        for i in 0..3 {
            let f = fs[i];
            let f_poly = Polynomial { coeffs: f.to_vec() };
            let poly_for_commit =
                compute_poly_coeff_round1(B_RANDS[2 * i], B_RANDS[2 * i + 1], f_poly);
            println!("{:?}", poly_for_commit);

            println!("Commited Point: {:?}", commit_poly(poly_for_commit));
        }
    }
}
