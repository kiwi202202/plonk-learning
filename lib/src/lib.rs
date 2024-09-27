pub mod field;
pub mod field_extension;
pub mod polynomial;
pub mod prescribed_permutation;
pub mod pythagorean_circuit;
pub mod pythagorean_transcript;
pub mod round1;
pub mod round2;
pub mod round3;
pub mod round4;
pub mod round5;
pub mod srs;
pub mod verifier;

pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

pub fn greet() {
    println!("Hello from the library!");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
