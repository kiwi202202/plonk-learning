use crate::field::F17;

#[derive(Debug, Clone)]
pub struct Polynomial {
    pub coeffs: Vec<F17>,
}

impl Polynomial {
    pub fn add(&self, other: &Polynomial) -> Polynomial {
        let max_len = usize::max(self.coeffs.len(), other.coeffs.len());
        let mut result = vec![F17::ZERO; max_len];

        for i in 0..self.coeffs.len() {
            result[i] = result[i].add(self.coeffs[i]);
        }

        for i in 0..other.coeffs.len() {
            result[i] = result[i].add(other.coeffs[i]);
        }

        Polynomial { coeffs: result }
    }

    // ax^b
    pub fn mul_by_monomial(&self, a: F17, b: usize) -> Polynomial {
        let mut result = vec![F17::ZERO; self.coeffs.len() + b];
        for i in 0..self.coeffs.len() {
            result[i + b] = self.coeffs[i].mul(a);
        }
        Polynomial { coeffs: result }
    }

    // Evaluate the polynomial at a given point `x`
    // Horner’s method
    pub fn evaluate(&self, x: F17) -> F17 {
        let mut result = F17::ZERO;
        // From the highest degree coeff, Step-down calculation
        for i in (0..self.coeffs.len()).rev() {
            result = result.mul(x).add(self.coeffs[i]);
        }
        result
    }

    // Multiply two polynomials
    pub fn mul(&self, other: &Polynomial) -> Polynomial {
        let mut result = vec![F17::ZERO; self.coeffs.len() + other.coeffs.len() - 1];

        for i in 0..self.coeffs.len() {
            for j in 0..other.coeffs.len() {
                result[i + j] = result[i + j].add(self.coeffs[i].mul(other.coeffs[j]));
            }
        }

        Polynomial { coeffs: result }
    }

    pub fn evaluate_at_omega_x(&self, omega: F17) -> Polynomial {
        let mut result = vec![F17::ZERO; self.coeffs.len()];

        for (i, &coeff) in self.coeffs.iter().enumerate() {
            let omega_i = omega.pow(i as u32);
            result[i] = coeff.mul(omega_i);
        }

        Polynomial { coeffs: result }
    }

    pub fn sub(&self, other: &Polynomial) -> Polynomial {
        let max_len = usize::max(self.coeffs.len(), other.coeffs.len());
        let mut result = vec![F17::ZERO; max_len];

        for i in 0..self.coeffs.len() {
            result[i] = result[i].add(self.coeffs[i]);
        }

        for i in 0..other.coeffs.len() {
            result[i] = result[i].sub(other.coeffs[i]);
        }

        Polynomial { coeffs: result }
    }

    // // Synthetic division by Z_H(x) = x^4 - 1
    // pub fn synthetic_div_zh(&self) -> (Polynomial, Polynomial) {
    //     let mut remainder = self.coeffs.clone(); // Initialize the remainder as t(x)
    //     let mut quotient = vec![F17::ZERO; self.coeffs.len() - 4]; // Initialize the quotient

    //     // Start from the highest degree term and proceed to divide
    //     for i in (4..self.coeffs.len()).rev() {
    //         // The current highest degree coefficient
    //         let current_coeff = remainder[i];

    //         // The quotient term is the current highest degree coefficient
    //         quotient[i - 4] = current_coeff;

    //         // Update the remainder: propagate the effect of (x^4 - 1) to the left
    //         remainder[i] = F17::ZERO; // Clear the highest degree term
    //         remainder[i - 4] = remainder[i - 4].sub(current_coeff); // Subtract (multiplied by -1)
    //     }

    //     // Return the quotient and remainder
    //     (
    //         Polynomial { coeffs: quotient },
    //         Polynomial {
    //             coeffs: remainder[0..4].to_vec(),
    //         }, // Remainder consists of the last 4 terms
    //     )
    // }

    // Perform polynomial long division, returning the quotient and remainder
    pub fn long_div(&self, divisor: &Polynomial) -> (Polynomial, Polynomial) {
        let mut remainder = self.coeffs.clone(); // The remainder starts as the dividend's coefficients
        let divisor_degree = divisor.coeffs.len() - 1; // The degree of the divisor (highest exponent)
        let mut quotient = vec![F17::ZERO; self.coeffs.len() - divisor.coeffs.len() + 1]; // Initialize quotient

        // Continue while the degree of the remainder is greater than or equal to the divisor's degree
        while remainder.len() >= divisor.coeffs.len() {
            // The leading term of the quotient is the leading coefficient of the remainder
            // divided by the leading coefficient of the divisor
            let lead_coeff_div = remainder
                .last()
                .unwrap()
                .mul(divisor.coeffs[divisor_degree].inv().unwrap());

            // The current term in the quotient has a degree difference between remainder and divisor
            let degree_diff = remainder.len() - divisor.coeffs.len();
            quotient[degree_diff] = lead_coeff_div;

            // Update the remainder by subtracting the result of (lead_coeff_div * divisor)
            for i in 0..divisor.coeffs.len() {
                let remainder_idx = degree_diff + i;
                remainder[remainder_idx] =
                    remainder[remainder_idx].sub(divisor.coeffs[i].mul(lead_coeff_div));
            }

            // Remove the highest-degree term from the remainder (it has been fully reduced)
            while remainder.last() == Some(&F17::ZERO) {
                remainder.pop();
            }
        }

        (
            Polynomial { coeffs: quotient },
            Polynomial { coeffs: remainder },
        )
    }

    pub fn split_into_three(&self) -> (Polynomial, Polynomial, Polynomial) {
        let total_len = self.coeffs.len();

        let low_len = total_len / 3;
        let mid_len = total_len / 3;
        let high_len = total_len - low_len - mid_len;

        let low_coeffs = self.coeffs[0..low_len].to_vec();

        let mid_coeffs = self.coeffs[low_len..(low_len + mid_len)].to_vec();

        let high_coeffs = self.coeffs[(low_len + mid_len)..total_len].to_vec();

        (
            Polynomial { coeffs: low_coeffs },
            Polynomial { coeffs: mid_coeffs },
            Polynomial {
                coeffs: high_coeffs,
            },
        )
    }
}

#[allow(non_snake_case)]
pub fn get_Z_H() -> Polynomial {
    Polynomial {
        coeffs: vec![F17::NEG_ONE, F17::ZERO, F17::ZERO, F17::ZERO, F17::ONE], // -1 + x^4
    }
}
// pub const Z_H: Polynomial = ;
