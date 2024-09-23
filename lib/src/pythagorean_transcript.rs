use crate::field::F101;

pub fn gen_transcript() -> ([F101; 4], [F101; 4], [F101; 4]) {
    let a: [F101; 4] = [F101::new(3), F101::new(4), F101::new(5), F101::new(9)];
    let b: [F101; 4] = [F101::new(3), F101::new(4), F101::new(5), F101::new(16)];
    let c: [F101; 4] = [F101::new(9), F101::new(16), F101::new(25), F101::new(25)];
    return (a, b, c);
}
