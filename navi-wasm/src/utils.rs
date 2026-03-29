
use wasm_bindgen::prelude::*;
use nalgebra::DMatrix;

#[wasm_bindgen]
extern "C" {
    pub fn alert(s: &str);
}

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = console)]
    pub fn log(s: &str);
}

#[macro_export]
macro_rules! console_log {
    ($($t:tt)*) => (utils::log(&format_args!($($t)*).to_string()))
}

pub fn set_panic_hook() {
    #[cfg(feature = "console_error_panic_hook")]
    console_error_panic_hook::set_once();
}

pub fn new_square_dmatrix(v: &[f64]) -> Result<DMatrix<f64>, String> {
    let n2 = v.len();
    let n = n2.isqrt();
    if n * n != n2 {
        return Err("Input slice length is not a perfect square".to_string());
    }
    let a : DMatrix<f64> = DMatrix::from_row_slice(n, n, v);
    Ok(a)
}

pub fn new_dmatrix(v: &js_sys::Array<js_sys::Array<JsValue>>) -> Result<DMatrix<f64>, String> {

    let n = v.length() as usize;
    if n==0 {
        return Ok(DMatrix::from_element(0, 0, 0.0));
    }

    let row0 = &v.get(0 as u32);
    let m = row0.length() as usize;

    let mut a : DMatrix<f64> = DMatrix::from_element(n, m, 0.0);
    for i in 0..n {
        let row = &v.get(i as u32);
        if row.length() as usize != m {
            return Err(format!("Row {} has different length than the first row", i));
        }
        for j in 0..m {
            let val = row.get_checked(j as u32)
                .and_then(|v| v.as_f64())
                .ok_or_else(|| format!("Element ({}, {}) is not a number", i, j))?;
            a[(i, j)] = val;
        }
    }
    Ok(a)
}

pub fn js_value_from_slice(m: &[f64]) -> JsValue {
    let js_res = js_sys::Array::new();
    for i in 0..m.len() {
        js_res.push(&JsValue::from_f64(m[i]));
    }
    js_res.into()
}

pub fn new_js_value(m: &DMatrix<f64>) -> JsValue {
    let js_res = js_sys::Array::new();
    for i in 0..m.nrows() {
        let js_row = js_sys::Array::new();
        for j in 0..m.ncols() {
            js_row.push(&JsValue::from_f64(m[(i, j)]));
        }
        js_res.push(&js_row);
    }
    js_res.into()
}