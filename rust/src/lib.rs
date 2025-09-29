use pyo3::prelude::*;

/// Compute pairwise squared distances between two point sets (m x 3) and (n x 3).
/// Returns a Python list of floats flattened row-major (m*n) for initial prototype simplicity.
#[pyfunction]
fn pairwise_sq_dists(a: Vec<f32>, b: Vec<f32>, m: usize, n: usize) -> PyResult<Vec<f32>> {
    if a.len() != m*3 || b.len() != n*3 { return Err(pyo3::exceptions::PyValueError::new_err("Shape mismatch")); }
    let mut out = vec![0f32; m*n];
    for i in 0..m { for j in 0..n { let ai = i*3; let bj = j*3; let dx = a[ai]-b[bj]; let dy = a[ai+1]-b[bj+1]; let dz = a[ai+2]-b[bj+2]; out[i*n + j] = dx*dx + dy*dy + dz*dz; }}
    Ok(out)
}

#[pymodule]
fn molbridge_geom(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pairwise_sq_dists, m)?)?;
    Ok(())
}
