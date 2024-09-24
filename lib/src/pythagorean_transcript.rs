use crate::field::F17;

pub fn gen_transcript() -> ([F17; 4], [F17; 4], [F17; 4]) {
    let a: [F17; 4] = [F17::new(3), F17::new(4), F17::new(5), F17::new(9)];
    let b: [F17; 4] = [F17::new(3), F17::new(4), F17::new(5), F17::new(16)];
    let c: [F17; 4] = [F17::new(9), F17::new(16), F17::new(25), F17::new(25)];
    return (a, b, c);
}
