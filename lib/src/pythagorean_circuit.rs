// q_L * a + q_R * b + q_O * c + q_M * ab + q_C = 0
// a: input_left
// b: input_right
// c: output
// q_C: constant, used for public values

// We use four gates to express a simplified pythoagorean equation: d^2 + e^2 = f^2
// Gate 0: 0*a0 + 0*b0 + (-1)*c0 + 1*a0b0 + 0 = 0
// Gate 1: 0*a1 + 0*b1 + (-1)*c1 + 1*a1b1 + 0 = 0
// Gate 2: 0*a2 + 0*b2 + (-1)*c2 + 1*a2b2 + 0 = 0
// Gate 3: 1*a3 + 1*b3 + (-1)*c3 + 0*a3b3 + 0 = 0

use crate::{
    field::{solve_coefficients, F17},
    polynomial::Polynomial,
};

// Selectors, the verifier would have |S|
pub const Q_L: [F17; 4] = [F17::ZERO, F17::ZERO, F17::ZERO, F17::ONE];
pub const Q_R: [F17; 4] = [F17::ZERO, F17::ZERO, F17::ZERO, F17::ONE];
pub const Q_O: [F17; 4] = [F17::NEG_ONE, F17::NEG_ONE, F17::NEG_ONE, F17::NEG_ONE]; // -1 in field 101
pub const Q_M: [F17; 4] = [F17::ONE, F17::ONE, F17::ONE, F17::ZERO];
pub const Q_C: [F17; 4] = [F17::ZERO, F17::ZERO, F17::ZERO, F17::ZERO];

pub fn gen_s_polys_point_value() -> [[(F17, F17); 4]; 5] {
    [Q_L, Q_R, Q_O, Q_M, Q_C].map(|q| {
        F17::H
            .iter()
            .zip(q.iter())
            .map(|(&x, &y)| (x, y))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap()
    })
}

pub fn gen_s_polys_coeff() -> [[F17; 4]; 5] {
    let polys_pv = gen_s_polys_point_value();
    polys_pv
        .into_iter()
        .map(|poly_pv| solve_coefficients(&poly_pv).unwrap())
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
}

pub fn gen_selector_polys() -> [Polynomial; 5] {
    let [q_l, q_r, q_o, q_m, q_c] = gen_s_polys_coeff();

    [q_l, q_r, q_o, q_m, q_c].map(|coeffs| Polynomial {
        coeffs: coeffs.to_vec(),
    })
}

#[cfg(test)]
mod tests {
    use crate::field::solve_coefficients;

    use super::*;

    #[test]
    fn gen_s_polys_pv_test() {
        let polys_pv = gen_s_polys_point_value();
        for poly_pv in polys_pv {
            println!("poly in point-value form: {:?}", poly_pv);
            println!("poly in coeff form: {:?}", solve_coefficients(&poly_pv));
        }
    }

    #[test]
    fn gen_s_polys_coeff_test() {
        let polys_coeff = gen_s_polys_coeff();
        for poly_coeff in polys_coeff {
            println!("poly in coeff form: {:?}", poly_coeff);
        }
    }
}
