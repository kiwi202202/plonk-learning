use crate::{
    field::{Point, F17},
    polynomial::Polynomial,
    prescribed_permutation::gen_sigma_polys,
    pythagorean_circuit::gen_selector_polys,
    round1::gen_round1_result,
    round2::gen_round2_result,
    round3::gen_round3_result,
    round4::{gen_round4_result, Round4Output, ZETA},
    srs::commit_poly,
};

pub const V: F17 = F17(12);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Plonk_Proof {
    pub a_box: Point,
    pub b_box: Point,
    pub c_box: Point,
    pub z_box: Point,
    pub t_low_box: Point,
    pub t_mid_box: Point,
    pub t_high_box: Point,
    pub w_zeta_box: Point,
    pub w_zeta_omega_box: Point,
    pub a_bar: F17,
    pub b_bar: F17,
    pub c_bar: F17,
    pub sigma1_bar: F17,
    pub sigma2_bar: F17,
    pub r_bar: F17,
    pub z_omega_bar: F17,
}

pub fn gen_round5_result() -> Plonk_Proof {
    let [q_l, q_r, q_o, q_m, q_c] = gen_selector_polys();
    let ([a, b, c], [a_box, b_box, c_box]) = gen_round1_result();
    let [sigma1, sigma2, sigma3] = gen_sigma_polys();
    let (z_x, z_box) = gen_round2_result();
    let (t, t_low, t_mid, t_high, t_low_box, t_mid_box, t_high_box) = gen_round3_result();
    // let round4_output = gen_round4_result();
    let Round4Output {
        a_bar,
        b_bar,
        c_bar,
        sigma1_bar,
        sigma2_bar,
        t_bar,
        z_omega_bar,
        r_bar,
        r,
    } = gen_round4_result();
    // number of gates
    let n = 4u32;
    println!("t_low:{:?}", t_low);
    println!("t_mid_zeta:{:?}", t_mid.mul_by_monomial(ZETA.pow(n + 2), 0));
    // println!("t_high:{:?}", t_high);
    println!(
        "t_high_zeta:{:?}",
        t_high.mul_by_monomial(ZETA.pow(2 * n + 4), 0)
    );
    let term1 = t_low
        .add(&t_mid.mul_by_monomial(ZETA.pow(n + 2), 0))
        .add(&t_high.mul_by_monomial(ZETA.pow(2 * n + 4), 0))
        .sub(&Polynomial {
            coeffs: vec![t_bar],
        });
    let term2 = r
        .sub(&Polynomial {
            coeffs: vec![r_bar],
        })
        .mul_by_monomial(V, 0);
    let term3 = a
        .sub(&Polynomial {
            coeffs: vec![a_bar],
        })
        .mul_by_monomial(V.pow(2), 0);
    let term4 = b
        .sub(&Polynomial {
            coeffs: vec![b_bar],
        })
        .mul_by_monomial(V.pow(3), 0);
    let term5 = c
        .sub(&Polynomial {
            coeffs: vec![c_bar],
        })
        .mul_by_monomial(V.pow(4), 0);
    let term6 = sigma1
        .sub(&Polynomial {
            coeffs: vec![sigma1_bar],
        })
        .mul_by_monomial(V.pow(5), 0);
    let term7 = sigma2
        .sub(&Polynomial {
            coeffs: vec![sigma2_bar],
        })
        .mul_by_monomial(V.pow(6), 0);
    let w_zeta_x_neg_zeta = term1
        .add(&term2)
        .add(&term3)
        .add(&term4)
        .add(&term5)
        .add(&term6)
        .add(&term7);
    println!("w_zeta_x_neg_zeta: {:?}", w_zeta_x_neg_zeta);
    let (w_zeta, _remainder) = w_zeta_x_neg_zeta.long_div(&Polynomial {
        coeffs: vec![ZETA.neg(), F17::ONE],
    });
    println!("w_zeta: {:?}", w_zeta);

    // println!("z_x: {:?}", z_x);
    // let w_zeta_omega_nomi = z_x.sub(&Polynomial {
    //     coeffs: vec![z_omega_bar.neg()],
    // });
    // println!("w_zeta_omega_nomi: {:?}", w_zeta_omega_nomi);
    // println!("z_omega_bar_neg: {:?}", z_omega_bar.neg());

    let (w_zeta_omega, _remainder) = z_x
        .sub(&Polynomial {
            coeffs: vec![z_omega_bar],
        })
        .long_div(&Polynomial {
            coeffs: vec![ZETA.mul(F17::H[1]).neg(), F17::ONE],
        });
    println!("w_zeta_omega: {:?}", w_zeta_omega);

    let w_zeta_box = commit_poly(w_zeta);
    let w_zeta_omega_box = commit_poly(w_zeta_omega);
    println!("w_zeta_box: {:?}", w_zeta_box);
    println!("w_zeta_omega_box: {:?}", w_zeta_omega_box);
    Plonk_Proof {
        a_box,
        b_box,
        c_box,
        z_box,
        t_low_box,
        t_mid_box,
        t_high_box,
        w_zeta_box,
        w_zeta_omega_box,
        a_bar,
        b_bar,
        c_bar,
        sigma1_bar,
        sigma2_bar,
        r_bar,
        z_omega_bar,
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn gen_round5_result_test() {
        let proof = gen_round5_result();
        println!("plonk proof: {:?}", proof);
    }
}
