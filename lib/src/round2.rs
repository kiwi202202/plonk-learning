use crate::{
    field::{solve_coefficients, Point, F17},
    polynomial::{get_Z_H, Polynomial},
    prescribed_permutation::gen_W_coeff,
    pythagorean_transcript::gen_transcript,
    srs::commit_poly,
};

// generate random b1..b9 in F17 (b1..b6 is from round1)
// Prover roll the dice and generate the random number
pub const B_RANDS: [F17; 9] = [
    F17(7),
    F17(4),
    F17(11),
    F17(12),
    F17(16),
    F17(2),
    F17(14),
    F17(11),
    F17(7),
];

// Challenges
// Verifier roll the dice and generate the random number
pub const BETA: F17 = F17(12);
pub const GAMMA: F17 = F17(13);

pub fn compute_acc(
    acc_prev: F17,
    i: usize,
    a: [F17; 4],
    b: [F17; 4],
    c: [F17; 4],
    s1: [F17; 4],
    s2: [F17; 4],
    s3: [F17; 4], // S_sigma
    omega: [F17; 4],
) -> F17 {
    let ai = a[i - 1];
    let bi = b[i - 1];
    let ci = c[i - 1];

    println!("ai,bi,ci: {:?}, {:?}, {:?}", ai, bi, ci);
    // ω^(i-1)
    let omega_i_neg_one = omega[i - 1];

    let s1_eval = Polynomial {
        coeffs: s1.to_vec(),
    }
    .evaluate(omega_i_neg_one);
    let s2_eval = Polynomial {
        coeffs: s2.to_vec(),
    }
    .evaluate(omega_i_neg_one);
    let s3_eval = Polynomial {
        coeffs: s3.to_vec(),
    }
    .evaluate(omega_i_neg_one);

    let k1 = F17(2);
    let k2 = F17(3);

    let numerator1 = ai.add(BETA.mul(omega_i_neg_one)).add(GAMMA);
    let numerator2 = bi.add(BETA.mul(k1).mul(omega_i_neg_one)).add(GAMMA);
    let numerator3 = ci.add(BETA.mul(k2).mul(omega_i_neg_one)).add(GAMMA);

    let denominator1 = ai.add(BETA.mul(s1_eval)).add(GAMMA);
    let denominator2 = bi.add(BETA.mul(s2_eval)).add(GAMMA);
    let denominator3 = ci.add(BETA.mul(s3_eval)).add(GAMMA);

    let numerator = numerator1.mul(numerator2).mul(numerator3);
    let denominator = denominator1.mul(denominator2).mul(denominator3);

    // acc_i = acc_(i-1) * (numerator / denominator)
    acc_prev.mul(numerator.mul(denominator.inv().unwrap()))
}

pub fn initial_acc() -> F17 {
    F17(1)
}

pub fn compute_poly_coeff_round2(b_xx: F17, b_x: F17, b: F17, f: Polynomial) -> Polynomial {
    let zh_poly = get_Z_H();
    // compute (b1 * x + b2) * Z_H(x)
    let after_bxx = zh_poly.mul_by_monomial(b_xx, 2); // b_xx * x^2 * Z_H(x)
    let after_bx = zh_poly.mul_by_monomial(b_x, 1); // b_x * x * Z_H(x)
    let after_b_const = zh_poly.mul_by_monomial(b, 0); // b * Z_H(x)
    let final_zh = after_bxx.add(&after_bx.add(&after_b_const)); // b1 * x * Z_H(x) + b2 * Z_H(x)

    // a(x) = (b1x + b2) * Z_H(x) + f_a(x)
    final_zh.add(&f)
}

pub fn gen_round2_result() -> (Polynomial, Point) {
    let acc0 = initial_acc();
    let (a, b, c) = gen_transcript();
    let (s1, s2, s3) = gen_W_coeff();
    let acc1 = compute_acc(acc0, 1, a, b, c, s1, s2, s3, F17::H);
    let acc2 = compute_acc(acc1, 2, a, b, c, s1, s2, s3, F17::H);
    let acc3 = compute_acc(acc2, 3, a, b, c, s1, s2, s3, F17::H);

    let acc_points = F17::H
        .into_iter()
        .zip([acc0, acc1, acc2, acc3].into_iter())
        .map(|(x, y)| (x, y))
        .collect::<Vec<_>>();

    let acc_coeff = solve_coefficients(&acc_points.try_into().unwrap());
    println!("acc coeff: {:?}", acc_coeff.unwrap());

    let round2_poly = compute_poly_coeff_round2(
        B_RANDS[6],
        B_RANDS[7],
        B_RANDS[8],
        Polynomial {
            coeffs: acc_coeff.unwrap().to_vec(),
        },
    );

    println!("round2 poly: {:?}", round2_poly);

    let commited_point = commit_poly(round2_poly.clone());
    println!("round2 committed point: {:?}", commited_point);
    (round2_poly, commited_point)
}

#[cfg(test)]
mod tests {
    use crate::{
        field::solve_coefficients, prescribed_permutation::gen_W_coeff,
        pythagorean_transcript::gen_transcript, srs::commit_poly,
    };

    use super::*;

    #[test]
    fn compute_acc_test() {
        let acc0 = initial_acc();
        println!("acc0: {:?}", acc0);
        let (a, b, c) = gen_transcript();
        let (s1, s2, s3) = gen_W_coeff();
        let acc1 = compute_acc(acc0, 1, a, b, c, s1, s2, s3, F17::H);
        println!("acc1: {:?}", acc1);
        let acc2 = compute_acc(acc1, 2, a, b, c, s1, s2, s3, F17::H);
        println!("acc2: {:?}", acc2);
        let acc3 = compute_acc(acc2, 3, a, b, c, s1, s2, s3, F17::H);
        println!("acc3: {:?}", acc3);

        let acc_points = F17::H
            .into_iter()
            .zip([acc0, acc1, acc2, acc3].into_iter())
            .map(|(x, y)| (x, y))
            .collect::<Vec<_>>();

        let acc_coeff = solve_coefficients(&acc_points.try_into().unwrap());
        println!("acc coeff: {:?}", acc_coeff.unwrap());

        let round2_poly = compute_poly_coeff_round2(
            B_RANDS[6],
            B_RANDS[7],
            B_RANDS[8],
            Polynomial {
                coeffs: acc_coeff.unwrap().to_vec(),
            },
        );

        println!("round2 poly: {:?}", round2_poly);

        let commited_point = commit_poly(round2_poly);
        println!("round2 committed point: {:?}", commited_point);
    }
}
