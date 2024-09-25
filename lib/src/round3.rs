use crate::{
    field::{Point, F17},
    polynomial::{get_Z_H, Polynomial},
    prescribed_permutation::gen_sigma_polys,
    pythagorean_circuit::{gen_s_polys_coeff, gen_selector_polys},
    round1::gen_round1_result,
    round2::{gen_round2_result, BETA, GAMMA},
    srs::commit_poly,
};

pub const ALPHA: F17 = F17(15);

// Langrange base, calculate by interpolating (1,0,0,0) on F17::H
pub const L1: [F17; 4] = [F17(13), F17(13), F17(13), F17(13)];

pub fn compute_poly_coeff_round3() -> (Polynomial, (Polynomial, Polynomial, Polynomial)) {
    let [q_l, q_r, q_o, q_m, q_c] = gen_selector_polys();
    let ([a, b, c], [a_box, b_box, c_box]) = gen_round1_result();
    let [sigma1, sigma2, sigma3] = gen_sigma_polys();
    let zh = get_Z_H();
    let a_b_qM = a.mul(&b).mul(&q_m);
    println!("a(x) * b(x) * qM(x): {:?}", a.mul(&b).mul(&q_m));

    let a_qL = a.mul(&q_l);
    println!("a(x) * qL(x): {:?}", a.mul(&q_l));

    let b_qR = b.mul(&q_r);
    println!("b(x) * qR(x): {:?}", b.mul(&q_r));

    let c_qO = c.mul(&q_o);
    println!("c(x) * qO(x): {:?}", c.mul(&q_o));

    let public_inputs_poly = Polynomial {
        coeffs: vec![F17::ZERO],
    };

    let beta_x = Polynomial {
        coeffs: [F17::ZERO, BETA].to_vec(),
    };

    let gamma_constant = Polynomial {
        coeffs: [GAMMA].to_vec(),
    };

    let alpha_a_beta_x_gamma = a
        .add(&beta_x)
        .add(&gamma_constant)
        .mul_by_monomial(ALPHA, 0);

    println!("alpha_a_beta_x_gamma: {:?}", alpha_a_beta_x_gamma);

    let k1 = F17(2);
    let k2 = F17(3);
    let b_beta_k1_x_gamma = b.add(&beta_x.mul_by_monomial(k1, 0)).add(&gamma_constant);
    println!("b_beta_k1_x_gamma: {:?}", b_beta_k1_x_gamma);

    let c_beta_k2_x_gamma = c.add(&beta_x.mul_by_monomial(k2, 0)).add(&gamma_constant);
    println!("c_beta_k2_x_gamma: {:?}", c_beta_k2_x_gamma);

    let (z_x, z_box) = gen_round2_result();
    println!("z_x: {:?}", z_x);

    let z_omega_x = z_x.evaluate_at_omega_x(F17::H[1]);
    println!("z_omega_x: {:?}", z_omega_x);

    let alpha_a_beta_sigma1_gamma = a
        .add(&sigma1.mul_by_monomial(BETA, 0))
        .add(&gamma_constant)
        .mul_by_monomial(ALPHA, 0);

    println!("alpha_a_beta_sigma1_gamma: {:?}", alpha_a_beta_sigma1_gamma);

    let b_beta_k1_sigma2_gamma = b.add(&sigma2.mul_by_monomial(BETA, 0)).add(&gamma_constant);
    println!("b_beta_k1_sigma2_gamma: {:?}", b_beta_k1_sigma2_gamma);

    let c_beta_k2_sigma3_gamma = c.add(&sigma3.mul_by_monomial(BETA, 0)).add(&gamma_constant);
    println!("c_beta_k2_sigma3_gamma: {:?}", c_beta_k2_sigma3_gamma);

    let z_x_neg_one = z_x.add(&Polynomial {
        coeffs: vec![F17::NEG_ONE],
    });

    let alpha_2_z_x_negone_L1 = z_x_neg_one
        .mul(&Polynomial {
            coeffs: L1.to_vec(),
        })
        .mul_by_monomial(ALPHA.mul(ALPHA), 0);

    println!("alpha_2_z_x_negone_L1: {:?}", alpha_2_z_x_negone_L1);

    let term1 = a_b_qM
        .add(&a_qL)
        .add(&b_qR)
        .add(&c_qO)
        .add(&public_inputs_poly)
        .add(&q_c);
    println!("term1: {:?}", term1);
    let term2 = alpha_a_beta_x_gamma
        .mul(&b_beta_k1_x_gamma)
        .mul(&c_beta_k2_x_gamma)
        .mul(&z_x);
    println!("term2: {:?}", term2);
    let term3 = alpha_a_beta_sigma1_gamma
        .mul(&b_beta_k1_sigma2_gamma)
        .mul(&c_beta_k2_sigma3_gamma)
        .mul(&z_omega_x);
    println!("term3: {:?}", term3);
    let term4 = alpha_2_z_x_negone_L1;
    println!("term4: {:?}", term4);

    let t_zh = term1.add(&term2).sub(&term3).add(&term4);
    println!("t_zh: {:?}", t_zh);

    let (t, t_remainder) = t_zh.long_div(&zh);
    println!("t: {:?}, t_remainder: {:?}", t, t_remainder);

    (t.clone(), t.split_into_three())
}

pub fn gen_round3_result() -> (
    Polynomial,
    Polynomial,
    Polynomial,
    Polynomial,
    Point,
    Point,
    Point,
) {
    let (t, (t_low, t_mid, t_high)) = compute_poly_coeff_round3();
    let t_low_box = commit_poly(t_low.clone());
    let t_mid_box = commit_poly(t_mid.clone());
    let t_high_box = commit_poly(t_high.clone());
    (t, t_low, t_mid, t_high, t_low_box, t_mid_box, t_high_box)
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn compute_poly_coeff_round3_test() {
        compute_poly_coeff_round3();
    }
}
