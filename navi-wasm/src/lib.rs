mod utils;

use wasm_bindgen::prelude::*;
use js_sys::Float64Array;

// https://www.nalgebra.rs/docs/user_guide/wasm_and_embedded_targets

#[wasm_bindgen]
pub fn greet() {
    utils::set_panic_hook();
    // utils::alert("Hello, navi-wasm!");
    console_log!("Hello, {}!", "navi-wasm");
}

#[wasm_bindgen]
pub fn vadd(a: &[f64]) -> f64 {
    utils::set_panic_hook();
    let mut s = 0.0;
    for x in a.iter() {
        s += x;
    }
    s
}

#[wasm_bindgen]
pub fn vadd2(a: &mut [f64]) {
    utils::set_panic_hook();
    for x in a.iter_mut() {
        *x = *x + 2.0;
    }
}

#[wasm_bindgen]
pub fn vadd3(a0: &[f64]) -> Float64Array {
    utils::set_panic_hook();
    let mut a = a0.to_vec();
    for x in a.iter_mut() {
        *x = *x + 2.0
    }
    Float64Array::from(&a[..])
}

#[wasm_bindgen]
pub fn cholesky(v: &[f64]) -> Result<JsValue, String>  {
    utils::set_panic_hook();
    let mtx = utils::new_square_dmatrix(v)?;
    let Some(decomp) = mtx.cholesky() else {
        return Err("Cholesky decomposition failed: matrix is not positive definite".to_string());
    };
    let res = decomp.l();
    // Ok(Float64Array::from(res.as_slice()))
    Ok(utils::new_js_value(&res))
}

#[wasm_bindgen]
pub fn svd(m: &js_sys::Array<js_sys::Array<JsValue>>) -> Result<JsValue, String>  {
    utils::set_panic_hook();
    let mtx= utils::new_dmatrix(m)?;
    let svd = mtx.svd(true, true);
    let Some(u) = svd.u else {
        return Err("SVD decomposition failed: could not compute U matrix".to_string());
    };
    let s = svd.singular_values;
    let Some(v_t) = svd.v_t else {
        return Err("SVD decomposition failed: could not compute V^T matrix".to_string());
    };

    let s2 = s.data.as_slice();

    // Map and object are both possible, but object is more convenient to use in JS

    // let res = js_sys::Map::new();
    // res.set(&"u".into(), &utils::new_js_value(&u));
    // res.set(&"s".into(), &utils::new_js_value(&s));
    // res.set(&"v_t".into(), &utils::new_js_value(&v_t));

    let res = js_sys::Object::new();
    let _ = js_sys::Reflect::set(&res, &"u".into(), &utils::new_js_value(&u));
    let _ = js_sys::Reflect::set(&res, &"s".into(), &utils::js_value_from_slice(s2))  ;
    let _ = js_sys::Reflect::set(&res, &"v_t".into(), &utils::new_js_value(&v_t));

    Ok(res.into())
}

#[derive(serde::Serialize, serde::Deserialize)]
struct NaviFLow {
    m: usize,
    n: usize,
    r: f64,
}

#[derive(serde::Serialize, serde::Deserialize)]
struct NaviFLowResult {
    grid: Vec<Vec<f64>>,
}

#[wasm_bindgen]
pub fn potential_flow(js_args: JsValue) -> Result<JsValue, String>  {
    let args = serde_wasm_bindgen::from_value::<NaviFLow>(js_args);
    let prms = match args {
        Ok(prms) => prms,
        Err(err) => return Err(format!("Failed to parse input: {}", err)),
    };

    let c = ( prms.n/2 ) as f64;
    let r2 = prms.r * prms.r;

    // let geom = navi_rs::pflow::FreeFlow {};
    let geom = navi_rs::pflow::Cylinder { center: nalgebra::Vector2::new(c, c), radius2: r2 };

    let grid = navi_rs::pflow::finite_diff(prms.n, prms.m, &geom);
    let data = utils::to_2d_vec(&grid);

    let res = NaviFLowResult { grid: data };
    let json = serde_wasm_bindgen::to_value(&res);
    match json {
        Ok(json) => Ok(json),
        Err(e) => Err(format!("Failed to serialize to JSON: {}", e)),
    }
}