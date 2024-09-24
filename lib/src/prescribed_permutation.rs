use crate::field::{solve_coefficients, F17};

// imagine a one-place rotation on all the wires
// left inputs are encoded by F17.H
// right inputs are encoded by F17.k1H
// outputs are encoded by F17.k2H
// We have such copy constraints:
// a1=b1, a2=b2, a3=b3, a4=c1
// b1=a1, b2=a2, b3=a3, b4=c2
// c1=a4, c2=b4, c3=c4, c4=c3
#[allow(non_snake_case)]
pub fn gen_W_coeff() -> ([F17; 4], [F17; 4], [F17; 4]) {
    let sigma_L = [F17(2), F17(8), F17(15), F17(3)];
    let sigma_R = [F17(1), F17(4), F17(16), F17(12)];
    let sigma_O = [F17(13), F17(9), F17(5), F17(14)];
    let sigma_L_points: [(F17, F17); 4] = F17::H
        .into_iter()
        .zip(sigma_L.into_iter())
        .map(|(x, y)| (x, y))
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    let sigma_R_points: [(F17, F17); 4] = F17::H
        .into_iter()
        .zip(sigma_R.into_iter())
        .map(|(x, y)| (x, y))
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    let sigma_O_points: [(F17, F17); 4] = F17::H
        .into_iter()
        .zip(sigma_O.into_iter())
        .map(|(x, y)| (x, y))
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    let sigma_L_coeff = solve_coefficients(&sigma_L_points).unwrap();
    let sigma_R_coeff = solve_coefficients(&sigma_R_points).unwrap();
    let sigma_O_coeff = solve_coefficients(&sigma_O_points).unwrap();
    (sigma_L_coeff, sigma_R_coeff, sigma_O_coeff)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gen_W_coeff_test() {
        let result = gen_W_coeff();
        println!("{:?}", result);
    }
}
