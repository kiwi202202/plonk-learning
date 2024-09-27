// The verifier has |W| (permutation commitment) and |S| (selector commitment) in advance.
// |W| is sigma1,2,3, |S| is q_L,q_R,q_O,q_M,q_C

use crate::{
    field::{point_add, scalar_mult, Point, F101, F17},
    prescribed_permutation::gen_sigma_polys,
    pythagorean_circuit::{gen_selector_polys, N},
    round2::{BETA, GAMMA},
    round3::ALPHA,
    round4::ZETA,
    round5::{Plonk_Proof, V},
    srs::commit_poly,
};

fn verify_point_on_ec_curve(p: Point) -> bool {
    let result = match p {
        Point::Infinity => true,
        Point::Point { x, y } => {
            let lhs = y.mul(y);
            let rhs = x.mul(x).mul(x).add(F101(3));
            lhs == rhs
        }
    };
    result
}

#[allow(non_snake_case)]
fn verify_ele_in_F17(ele: F17) -> bool {
    let num = ele.0;
    return num < 17u32;
}

pub fn verifier_process(proof: Plonk_Proof) {
    let [q_l, q_r, q_o, q_m, q_c] = gen_selector_polys();
    let q_l_box = commit_poly(q_l);
    let q_r_box = commit_poly(q_r);
    let q_o_box = commit_poly(q_o);
    let q_m_box = commit_poly(q_m);
    let q_c_box = commit_poly(q_c);

    println!(
        "{:?}, {:?}, {:?}, {:?}, {:?}",
        q_l_box, q_r_box, q_o_box, q_m_box, q_c_box
    );

    let [sigma1, sigma2, sigma3] = gen_sigma_polys();
    let sigma1_box = commit_poly(sigma1);
    let sigma2_box = commit_poly(sigma2);
    let sigma3_box = commit_poly(sigma3);
    println!("{:?}, {:?},  {:?}", sigma1_box, sigma2_box, sigma3_box);

    let random_u = F17(4);

    let Plonk_Proof {
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
    } = proof;

    // Step 1: check all the commitments are valid ellptic curve elements
    assert!(verify_point_on_ec_curve(a_box));
    assert!(verify_point_on_ec_curve(b_box));
    assert!(verify_point_on_ec_curve(c_box));
    assert!(verify_point_on_ec_curve(z_box));
    assert!(verify_point_on_ec_curve(t_low_box));
    assert!(verify_point_on_ec_curve(t_mid_box));
    assert!(verify_point_on_ec_curve(t_high_box));
    assert!(verify_point_on_ec_curve(w_zeta_box));
    assert!(verify_point_on_ec_curve(w_zeta_omega_box));

    // Step 2: check all the commitments are valid ellptic curve elements
    assert!(verify_ele_in_F17(a_bar));
    assert!(verify_ele_in_F17(b_bar));
    assert!(verify_ele_in_F17(c_bar));
    assert!(verify_ele_in_F17(sigma1_bar));
    assert!(verify_ele_in_F17(sigma2_bar));
    assert!(verify_ele_in_F17(r_bar));
    assert!(verify_ele_in_F17(z_omega_bar));

    // Step 3: check w_{i /in public input set} is valid F17 elements
    // skip, we have no public inputs for now

    // Step 4: compute Zeta^n - 1, n is the # of gates, 4
    // log(n) computation here
    // let n = 4u32;
    let z_h_eval = ZETA.pow(N).sub(F17::ONE);
    println!("z_h_eval:{:?}", z_h_eval);

    // Step 5:
    let l_1_eval = z_h_eval.mul(F17(N).mul(ZETA.sub(F17::ONE)).inv().unwrap());
    println!("l_1_eval:{:?}", l_1_eval);

    // Step 6:
    // skip for now, we have no public inputs

    // Step 7: compute quotient polynomial evaluation
    let a_bar_beta_sigma1_gamma = a_bar.add(BETA.mul(sigma1_bar)).add(GAMMA);
    let b_bar_beta_sigma2_gamma = b_bar.add(BETA.mul(sigma2_bar)).add(GAMMA);
    let c_bar_gamma_z_omega_alpha = c_bar.add(GAMMA).mul(z_omega_bar).mul(ALPHA);
    let t_bar = r_bar
        .sub(
            a_bar_beta_sigma1_gamma
                .mul(b_bar_beta_sigma2_gamma)
                .mul(c_bar_gamma_z_omega_alpha),
        )
        .sub(l_1_eval.mul(ALPHA).mul(ALPHA))
        .mul(z_h_eval.inv().unwrap());
    println!("t_bar:{:?}", t_bar);

    // Step 8: the first part of batch polynomial commitment
    let mut term1 = scalar_mult(a_bar.mul(b_bar).mul(V).0, q_m_box);
    term1 = point_add(term1, scalar_mult(a_bar.mul(V).0, q_l_box));
    term1 = point_add(term1, scalar_mult(b_bar.mul(V).0, q_r_box));
    term1 = point_add(term1, scalar_mult(c_bar.mul(V).0, q_o_box));
    term1 = point_add(term1, scalar_mult(V.0, q_c_box));

    let term2 = scalar_mult(
        a_bar
            .add(BETA.mul(ZETA))
            .add(GAMMA)
            .mul(b_bar.add(BETA.mul(ZETA).mul(F17::K1)).add(GAMMA))
            .mul(c_bar.add(BETA.mul(ZETA).mul(F17::K2)).add(GAMMA))
            .mul(ALPHA)
            .mul(V)
            .add(l_1_eval.mul(ALPHA).mul(ALPHA).mul(V))
            .add(random_u)
            .0,
        z_box,
    );

    let term3 = scalar_mult(
        a_bar
            .add(BETA.mul(sigma1_bar))
            .add(GAMMA)
            .mul(b_bar.add(BETA.mul(sigma2_bar)).add(GAMMA))
            .mul(ALPHA)
            .mul(V)
            .mul(BETA)
            .mul(z_omega_bar)
            .0,
        sigma3_box,
    )
    .point_neg();

    let d_box = term1.add(term2).add(term3);
    println!("d_box:{:?}", d_box);

    // q_m_box. a_bar
    //     .mul(b_bar)
    //     .mul()
    //     .mul(V)
    //     .add(a_bar.mul(V).mul(q_l_box))
    //     .add(b_bar.mul(V).mul(q_r_box))
    //     .add(c_bar.mul(V).mul(q_o_box))
    //     .add(V.mul(q_c_box));

    // Step 9: compute full batched polynomial commitment
    let f_box = t_low_box
        .add(t_mid_box.mul(ZETA.pow(N + 2)))
        .add(t_high_box.mul(ZETA.pow(2 * N + 4)))
        .add(d_box)
        .add(a_box.mul(V.pow(2)))
        .add(b_box.mul(V.pow(3)))
        .add(c_box.mul(V.pow(4)))
        .add(sigma1_box.mul(V.pow(5)))
        .add(sigma2_box.mul(V.pow(6)));
    println!("f_box:{:?}", f_box);

    // Step 10: compute group encoded batch evaluation
    let point_primitive = Point::Point {
        x: F101(1),
        y: F101(2),
    };
    let e_coeff = t_bar
        .add(V.mul(r_bar))
        .add(V.pow(2).add(a_bar))
        .add(V.pow(3).add(b_bar))
        .add(V.pow(4).add(c_bar))
        .add(V.pow(5).add(sigma1_bar))
        .add(V.pow(6).add(sigma2_bar))
        .add(random_u.mul(z_omega_bar));
    let e_box = point_primitive.mul(e_coeff);
    println!("e_box:{:?}", e_box);

    // Step 11: final pairing
}

#[cfg(test)]
mod tests {
    use crate::round5::gen_round5_result;

    use super::*;

    #[test]
    fn verifier_process_test() {
        let proof = gen_round5_result();
        verifier_process(proof);
    }
}
