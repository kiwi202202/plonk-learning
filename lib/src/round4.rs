use crate::{
    field::F17,
    polynomial::{get_Z_H, Polynomial},
    prescribed_permutation::gen_sigma_polys,
    pythagorean_circuit::gen_selector_polys,
    round1::gen_round1_result,
    round2::{gen_round2_result, BETA, GAMMA},
    round3::{gen_round3_result, ALPHA, L1},
};

pub const ZETA: F17 = F17(5);

pub struct Round4Output {
    pub a_bar: F17,
    pub b_bar: F17,
    pub c_bar: F17,
    pub sigma1_bar: F17,
    pub sigma2_bar: F17,
    pub t_bar: F17,
    pub z_omega_bar: F17,
    pub r_bar: F17,
    pub r: Polynomial,
}

pub fn gen_round4_result() -> Round4Output {
    let [q_l, q_r, q_o, q_m, q_c] = gen_selector_polys();
    let ([a, b, c], [a_box, b_box, c_box]) = gen_round1_result();
    let [sigma1, sigma2, sigma3] = gen_sigma_polys();
    let (t, t_low, t_mid, t_high, t_low_box, t_mid_box, t_high_box) = gen_round3_result();
    let (z_x, z_box) = gen_round2_result();
    println!("z_x: {:?}", z_x);

    let z_omega_x = z_x.evaluate_at_omega_x(F17::H[1]);
    println!("z_omega_x: {:?}", z_omega_x);

    let zh = get_Z_H();
    let a_bar = a.evaluate(ZETA);
    let b_bar = b.evaluate(ZETA);
    let c_bar = c.evaluate(ZETA);
    println!("a_bar: {:?}", a_bar);
    println!("b_bar: {:?}", b_bar);
    println!("c_bar: {:?}", c_bar);
    let sigma1_bar = sigma1.evaluate(ZETA);
    let sigma2_bar = sigma2.evaluate(ZETA);
    // let sigma3_bar = sigma3.evaluate(ZETA);
    println!("sigma1_bar: {:?}", sigma1_bar);
    println!("sigma2_bar: {:?}", sigma2_bar);
    // println!("sigma3_bar: {:?}", sigma3_bar);

    let t_bar = t.evaluate(ZETA);
    println!("t_bar: {:?}", t_bar);

    let z_omega_bar = z_omega_x.evaluate(ZETA);
    println!("z_omega_bar: {:?}", z_omega_bar);

    // term 1
    let ab_bar_q_m = q_m.mul_by_monomial(a_bar, 0).mul_by_monomial(b_bar, 0);
    println!("ab_bar_q_m: {:?}", ab_bar_q_m);

    let a_bar_q_l = q_l.mul_by_monomial(a_bar, 0);
    let b_bar_q_r = q_r.mul_by_monomial(b_bar, 0);
    let c_bar_q_o = q_o.mul_by_monomial(c_bar, 0);
    println!("a_bar_q_l: {:?}", a_bar_q_l);
    println!("b_bar_q_r: {:?}", b_bar_q_r);
    println!("c_bar_q_o: {:?}", c_bar_q_o);
    println!("q_c: {:?}", q_c);

    let a_bar_beta_zeta_gamma = a_bar.add(BETA.mul(ZETA)).add(GAMMA);
    let b_bar_beta_k1_zeta_gamma = b_bar.add(BETA.mul(F17::K1).mul(ZETA)).add(GAMMA);
    let c_bar_beta_k2_zeta_gamma = c_bar.add(BETA.mul(F17::K2).mul(ZETA)).add(GAMMA);

    let term2 = z_x.mul_by_monomial(
        a_bar_beta_zeta_gamma
            .mul(b_bar_beta_k1_zeta_gamma)
            .mul(c_bar_beta_k2_zeta_gamma)
            .mul(ALPHA),
        0,
    );
    println!("term2: {:?}", term2);

    let a_bar_beta_sigma1_bar_gamma = a_bar.add(BETA.mul(sigma1_bar)).add(GAMMA);
    let b_bar_beta_sigma2_bar_gamma = b_bar.add(BETA.mul(sigma2_bar)).add(GAMMA);
    let term3 = sigma3.mul_by_monomial(
        a_bar_beta_sigma1_bar_gamma
            .mul(b_bar_beta_sigma2_bar_gamma)
            .mul(BETA)
            .mul(z_omega_bar)
            .mul(ALPHA),
        0,
    );
    println!("term3: {:?}", term3);

    let term4 = z_x.mul_by_monomial(
        Polynomial {
            coeffs: L1.to_vec(),
        }
        .evaluate(ZETA)
        .mul(ALPHA)
        .mul(ALPHA),
        0,
    );
    println!("term4: {:?}", term4);

    let r = ab_bar_q_m
        .add(&a_bar_q_l)
        .add(&b_bar_q_r)
        .add(&c_bar_q_o)
        .add(&q_c)
        .add(&term2)
        .add(&term3)
        .add(&term4);
    println!("r: {:?}", r);

    let r_bar = r.evaluate(ZETA);
    println!("r_bar: {:?}", r_bar);
    Round4Output {
        a_bar,
        b_bar,
        c_bar,
        sigma1_bar,
        sigma2_bar,
        t_bar,
        z_omega_bar,
        r_bar,
        r,
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn gen_round4_result_test() {
        gen_round4_result();
    }
}
