
fn main() {
    let expected = [
            ( 0.5, 38.2924922548026 ),
            ( 1.0, 68.2689492137086 ),
            ( 2.0, 95.4499736103642 ),
            ( 3.0, 99.7300203936740 ),
            ];

    for xp in expected {
        let x = xp.0;
        let p = 1.0 - 2.0*recipes::fcts::gaussian_cumulative(-x);
        let q = xp.1;
        println!("Gaussian cumulative with {:.2} stddev = {:.2}%. Expected = {:.2}%. Error = {:+e}%.  ", x, p*100.0, q, (p * 100.0 - q).abs());
    }
}
