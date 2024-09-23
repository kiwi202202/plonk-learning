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

// use crate::field::F101::ZERO;

use crate::field::F101;

// Selectors, the verifier would have |S|
pub const Q_L: [F101; 4] = [F101::ZERO, F101::ZERO, F101::ZERO, F101::ONE];
pub const Q_R: [F101; 4] = [F101::ZERO, F101::ZERO, F101::ZERO, F101::ONE];
pub const Q_O: [F101; 4] = [F101::NEG_ONE, F101::NEG_ONE, F101::NEG_ONE, F101::NEG_ONE]; // -1 in field 101
pub const Q_M: [F101; 4] = [F101::ONE, F101::ONE, F101::ONE, F101::ZERO];
pub const Q_C: [F101; 4] = [F101::ZERO, F101::ZERO, F101::ZERO, F101::ZERO];
