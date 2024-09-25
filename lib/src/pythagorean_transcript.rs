use crate::field::{solve_coefficients, F17};

pub fn gen_transcript() -> ([F17; 4], [F17; 4], [F17; 4]) {
    let a: [F17; 4] = [F17::new(3), F17::new(4), F17::new(5), F17::new(9)];
    let b: [F17; 4] = [F17::new(3), F17::new(4), F17::new(5), F17::new(16)];
    let c: [F17; 4] = [F17::new(9), F17::new(16), F17::new(25), F17::new(25)];
    return (a, b, c);
}

pub fn gen_t_polys_point_value() -> [[(F17, F17); 4]; 3] {
    let (a, b, c) = gen_transcript();
    [a, b, c].map(|q| {
        F17::H
            .iter()
            .zip(q.iter())
            .map(|(&x, &y)| (x, y))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap()
    })
}

pub fn gen_t_polys_coeff() -> [[F17; 4]; 3] {
    let polys_pv = gen_t_polys_point_value();
    polys_pv
        .into_iter()
        .map(|poly_pv| solve_coefficients(&poly_pv).unwrap())
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
}

#[cfg(test)]
mod tests {
    use crate::field::solve_coefficients;

    use super::*;

    #[test]
    fn gen_t_polys_pv_test() {
        let polys_pv = gen_t_polys_point_value();
        for poly_pv in polys_pv {
            println!("poly in point-value form: {:?}", poly_pv);
            println!("poly in coeff form: {:?}", solve_coefficients(&poly_pv));
        }
    }

    #[test]
    fn gen_t_polys_coeff_test() {
        let polys_coeff = gen_t_polys_coeff();
        for poly_coeff in polys_coeff {
            println!("poly in coeff form: {:?}", poly_coeff);
        }
    }
}
