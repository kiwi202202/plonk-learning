use crate::field::F101;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct F1012 {
    real: F101,
    imag: F101,
}

impl F1012 {
    pub fn new(real: u32, imag: u32) -> Self {
        F1012 {
            real: F101::new(real),
            imag: F101::new(imag),
        }
    }
}
