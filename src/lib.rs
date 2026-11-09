//! Floating-point port of Pazit Haim-Kislev's exhaustive MATLAB algorithm.
//! Coordinates are (q_1,...,q_n,p_1,...,p_n); all indices are zero based.
use itertools::Itertools;
use nalgebra::{DMatrix as Matrix, DVector as Vector};
use serde::{Deserialize, Serialize};
use std::io::Write;
use std::process::{Command, Stdio};

pub type Vertices = Vec<Vec<f64>>;
#[derive(Debug, Deserialize, Serialize)]
pub struct Triangulation {
    pub delaunay: Vec<Vec<usize>>,
    pub hull: Vec<Vec<usize>>,
}
#[derive(Debug, Serialize, Deserialize)]
pub struct Geometry {
    pub vertices: Vertices,
    pub normals: Vertices,
    pub faces: Vec<Vec<usize>>,
    pub vertex_to_face: Vec<Vec<bool>>,
    pub heights: Vec<f64>,
    pub centroid: Vec<f64>,
}
#[derive(Debug, Serialize)]
pub struct ResultData {
    pub capacity: f64,
    pub minimum: f64,
    pub order: Vec<usize>,
    pub beta: Vec<f64>,
    pub supports: u64,
    pub orders: u64,
    pub singular_kkt: u64,
    pub status: &'static str,
    pub maximum: f64,
    pub maximizing_order: Vec<usize>,
    pub maximizing_beta: Vec<f64>,
    pub closure_residual: f64,
    pub normalization_residual: f64,
}

pub fn validate_vertices(k: &Vertices) -> Result<usize, String> {
    let d = k.first().map_or(0, Vec::len);
    if d == 0
        || d % 2 != 0
        || k.len() < d + 1
        || k.iter()
            .any(|r| r.len() != d || r.iter().any(|v| !v.is_finite()))
    {
        return Err(
            "expected at least d+1 finite vertices with equal, positive even dimension d".into(),
        );
    }
    let differences = Matrix::from_fn(k.len() - 1, d, |i, j| k[i + 1][j] - k[0][j]);
    let scale = differences.amax();
    if scale == 0. || !scale.is_finite() {
        return Err("zero or nonfinite vertex extent".into());
    }
    if (differences / scale)
        .svd(false, false)
        .singular_values
        .iter()
        .filter(|&&s| s > 1e-10)
        .count()
        < d
    {
        return Err("vertices are not full dimensional at tolerance 1e-10".into());
    }
    Ok(d)
}

/// Only Qhull's simplex indices cross the Python boundary.
pub fn triangulate(k: &Vertices) -> Result<Triangulation, String> {
    validate_vertices(k)?;
    let python = std::env::var("EHZ_PYTHON").unwrap_or_else(|_| "python3".into());
    let adapter = std::env::var("EHZ_QHULL_ADAPTER")
        .unwrap_or_else(|_| format!("{}/scripts/qhull_adapter.py", env!("CARGO_MANIFEST_DIR")));
    let mut child = Command::new(python)
        .arg(adapter)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| e.to_string())?;
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&serde_json::to_vec(k).map_err(|e| e.to_string())?)
        .map_err(|e| e.to_string())?;
    let output = child.wait_with_output().map_err(|e| e.to_string())?;
    if !output.status.success() {
        return Err(format!(
            "Qhull adapter: {}",
            String::from_utf8_lossy(&output.stderr)
        ));
    }
    serde_json::from_slice(&output.stdout).map_err(|e| e.to_string())
}
fn differences(k: &Vertices, face: &[usize]) -> Matrix<f64> {
    Matrix::from_fn(face.len() - 1, k[0].len(), |i, j| {
        k[face[i + 1]][j] - k[face[0]][j]
    })
}
/// Volume-weighted Delaunay simplex centroid, as in barycenter.m.
pub fn barycenter(k: &Vertices, simplices: &[Vec<usize>]) -> Result<Vec<f64>, String> {
    let d = validate_vertices(k)?;
    let mut c = vec![0.; d];
    let mut total = 0.;
    for s in simplices {
        if s.len() != d + 1 || s.iter().any(|&i| i >= k.len()) {
            return Err("invalid Delaunay simplex".into());
        }
        let volume = differences(k, s).determinant().abs();
        total += volume;
        for (j, coordinate) in c.iter_mut().enumerate() {
            *coordinate += volume * s.iter().map(|&i| k[i][j]).sum::<f64>();
        }
    }
    if total <= 0. || !total.is_finite() {
        return Err("zero/nonfinite triangulated volume".into());
    }
    for x in &mut c {
        *x /= total * (d + 1) as f64;
    }
    Ok(c)
}
/// MATLAB KtimesT row order: K varies fastest.
pub fn ktimes_t(k: &Vertices, t: &Vertices) -> Vertices {
    t.iter()
        .flat_map(|p| k.iter().map(move |q| q.iter().chain(p).copied().collect()))
        .collect()
}
pub fn calc_heights(k: &Vertices, faces: &[Vec<usize>], normals: &Vertices) -> Vec<f64> {
    faces
        .iter()
        .zip(normals)
        .map(|(f, n)| k[f[0]].iter().zip(n).map(|(x, y)| x * y).sum())
        .collect()
}
pub fn calc_normals_and_faces(k: &Vertices, hull: &[Vec<usize>]) -> Result<Geometry, String> {
    let d = validate_vertices(k)?;
    if hull
        .iter()
        .any(|f| f.len() != d || f.iter().any(|&i| i >= k.len()))
    {
        return Err("invalid hull simplex".into());
    }
    let used: Vec<_> = hull.iter().flatten().copied().sorted().dedup().collect();
    let vertices: Vertices = used.iter().map(|&i| k[i].clone()).collect();
    let mut simplices: Vec<Vec<usize>> = hull
        .iter()
        .map(|f| {
            f.iter()
                .map(|i| used.binary_search(i).unwrap())
                .sorted()
                .collect()
        })
        .collect();
    simplices.sort();
    let mut normals: Vertices = vec![];
    let mut faces = vec![];
    let mut memberships = vec![];
    for f in simplices {
        let a = differences(&vertices, &f);
        if a.clone()
            .svd(false, false)
            .singular_values
            .iter()
            .filter(|&&s| s > 1e-10)
            .count()
            < d - 1
        {
            continue;
        }
        let mut normal: Vec<f64> = (0..d)
            .map(|j| {
                let minor =
                    Matrix::from_fn(d - 1, d - 1, |r, c| a[(r, if c < j { c } else { c + 1 })]);
                (if j % 2 == 0 { -1. } else { 1. }) * minor.determinant()
            })
            .collect();
        let length = normal.iter().map(|v| v * v).sum::<f64>().sqrt();
        for x in &mut normal {
            *x /= length;
        }
        // Equivalent orientation test: find the first vertex outside the face plane.
        for v in &vertices {
            let side: f64 = v
                .iter()
                .zip(&vertices[f[0]])
                .zip(&normal)
                .map(|((x, y), n)| (x - y) * n)
                .sum();
            if side.abs() > 1e-10 {
                if side > 0. {
                    for x in &mut normal {
                        *x = -*x;
                    }
                }
                break;
            }
        }
        let index = normals
            .iter()
            .position(|n| n.iter().zip(&normal).all(|(a, b)| (a - b).abs() <= 1e-8));
        let index = index.unwrap_or_else(|| {
            normals.push(normal);
            faces.push(f.clone());
            normals.len() - 1
        });
        memberships.push((f, index));
    }
    if normals.is_empty() {
        return Err("no nondegenerate facets".into());
    }
    let mut vertex_to_face = vec![vec![false; normals.len()]; vertices.len()];
    for (f, j) in memberships {
        for i in f {
            vertex_to_face[i][j] = true;
        }
    }
    let heights = calc_heights(&vertices, &faces, &normals);
    Ok(Geometry {
        vertices,
        normals,
        faces,
        vertex_to_face,
        heights,
        centroid: vec![],
    })
}
pub fn geometry(k: &Vertices) -> Result<Geometry, String> {
    validate_vertices(k)?;
    // Translate/scale before Qhull and cofactor arithmetic; retain the volume centroid.
    let extent = k
        .iter()
        .flat_map(|r| r.iter().zip(&k[0]).map(|(x, y)| (x - y).abs()))
        .fold(0_f64, f64::max);
    let scaled: Vertices = k
        .iter()
        .map(|r| r.iter().zip(&k[0]).map(|(x, y)| (x - y) / extent).collect())
        .collect();
    let tri = triangulate(&scaled)?;
    let c = barycenter(&scaled, &tri.delaunay)?;
    let centered: Vertices = scaled
        .iter()
        .map(|r| r.iter().zip(&c).map(|(x, c)| x - c).collect())
        .collect();
    let tri = triangulate(&centered)?;
    let mut g = calc_normals_and_faces(&centered, &tri.hull)?;
    // A joggled simplex must still define a supporting plane of the original input.
    for (n, h) in g.normals.iter().zip(&g.heights) {
        if !h.is_finite()
            || *h <= 0.
            || centered
                .iter()
                .any(|v| v.iter().zip(n).map(|(x, y)| x * y).sum::<f64>() - h > RESIDUAL_TOL)
        {
            return Err("Qhull reconstruction is not a supporting facet with an interior origin at numerical tolerance".into());
        }
    }
    for v in &mut g.vertices {
        for x in v {
            *x *= extent;
        }
    }
    for h in &mut g.heights {
        *h *= extent;
    }
    g.centroid = c.iter().zip(&k[0]).map(|(x, y)| y + extent * x).collect();
    Ok(g)
}
pub fn omega(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() / 2;
    (0..n).map(|i| a[i] * b[n + i] - a[n + i] * b[i]).sum()
}
const RANK_TOL: f64 = 1e-12;
const RESIDUAL_TOL: f64 = 1e-9;

/// Row selection with pivoting and reorthogonalization, on scaled rows.
fn independent_rows(a: &Matrix<f64>) -> Vec<usize> {
    let mut columns: Vec<Vector<f64>> = (0..a.nrows())
        .map(|i| a.row(i).transpose().into_owned())
        .collect();
    let mut indices: Vec<usize> = (0..columns.len()).collect();
    let mut chosen = vec![];
    for j in 0..columns.len().min(a.ncols()) {
        let p = (j..columns.len())
            .max_by(|&i, &k| {
                columns[i]
                    .norm_squared()
                    .total_cmp(&columns[k].norm_squared())
            })
            .unwrap();
        columns.swap(j, p);
        indices.swap(j, p);
        let norm = columns[j].norm();
        if norm <= RANK_TOL {
            break;
        }
        chosen.push(indices[j]);
        let q = &columns[j] / norm;
        for c in columns.iter_mut().skip(j + 1) {
            // A second pass limits loss of orthogonality.
            for _ in 0..2 {
                let dot = q.dot(c);
                *c -= &q * dot;
            }
        }
    }
    chosen
}
#[derive(Debug)]
pub enum EqualitySystem {
    Inconsistent,
    Unique(Matrix<f64>, Vector<f64>),
    Multiple(Matrix<f64>, Vector<f64>),
}
pub fn prepare_aeq_and_beq(mut a: Matrix<f64>, mut b: Vector<f64>) -> EqualitySystem {
    for i in 0..a.nrows() {
        let scale = a.row(i).norm().max(b[i].abs());
        if scale > 0. {
            a.row_mut(i).scale_mut(1. / scale);
            b[i] /= scale;
        }
    }
    let rows = independent_rows(&a);
    let augmented = Matrix::from_fn(a.nrows(), a.ncols() + 1, |i, j| {
        if j == a.ncols() {
            b[i]
        } else {
            a[(i, j)]
        }
    });
    if rows.len() != independent_rows(&augmented).len() {
        return EqualitySystem::Inconsistent;
    }
    let reduced = Matrix::from_fn(rows.len(), a.ncols(), |i, j| a[(rows[i], j)]);
    let rhs = Vector::from_iterator(rows.len(), rows.iter().map(|&i| b[i]));
    if rows.len() == a.ncols() {
        EqualitySystem::Unique(reduced, rhs)
    } else {
        EqualitySystem::Multiple(reduced, rhs)
    }
}

/// Relative rank test and residual check, with no absolute determinant cutoff.
fn solve_nonsingular(a: &Matrix<f64>, b: &Vector<f64>) -> Option<Vector<f64>> {
    let svd = a.clone().svd(true, true);
    let largest = svd.singular_values.max();
    if !largest.is_finite() || largest == 0. || svd.singular_values.min() <= RANK_TOL * largest {
        return None;
    }
    let x = svd.solve(b, RANK_TOL * largest).ok()?;
    let residual = (a * &x - b).norm();
    if x.iter().any(|v| !v.is_finite())
        || residual > RESIDUAL_TOL * (1. + a.norm() * x.norm() + b.norm())
    {
        None
    } else {
        Some(x)
    }
}
fn feasible(a: &Matrix<f64>, b: &Vector<f64>, x: &Vector<f64>) -> bool {
    x.iter().all(|v| v.is_finite() && *v >= 0.)
        && (a * x - b).norm() <= RESIDUAL_TOL * (1. + a.norm() * x.norm() + b.norm())
}
fn nonnegative(mut x: Vector<f64>) -> Option<Vector<f64>> {
    if x.iter().any(|v| !v.is_finite() || *v < -RANK_TOL) {
        return None;
    }
    // Remove only roundoff-level negatives, then recheck the full constraints.
    for v in x.iter_mut() {
        *v = v.max(0.);
    }
    Some(x)
}
/// Positive spanning of normals is equivalent to boundedness of H-input with h>0.
/// Vertices of {t>=0: P^T t=0, sum(t)=1} suffice to check strict positive spanning.
fn validate_bounded(p: &Vertices) -> Result<(), String> {
    let m = p.len();
    let d = p[0].len();
    let a = Matrix::from_fn(d + 1, m, |i, j| if i == d { 1. } else { p[j][i] });
    let b = Vector::from_fn(d + 1, |i, _| if i == d { 1. } else { 0. });
    let span = Matrix::from_fn(d, m, |i, j| p[j][i]);
    let mut row_scaled = span.clone();
    for i in 0..d {
        let scale = row_scaled.row(i).norm();
        if scale > 0. {
            row_scaled.row_mut(i).scale_mut(1. / scale);
        }
    }
    if independent_rows(&row_scaled).len() != d {
        return Err("facet normals do not span the ambient space (unbounded input)".into());
    }
    let mut covered = vec![false; m];
    for size in 2..=(d + 1).min(m) {
        for support in (0..m).combinations(size) {
            let sub = Matrix::from_fn(d + 1, size, |i, j| a[(i, support[j])]);
            if let EqualitySystem::Unique(reduced, rhs) =
                prepare_aeq_and_beq(sub.clone(), b.clone())
            {
                if let Some(x) = solve_nonsingular(&reduced, &rhs).and_then(nonnegative) {
                    if feasible(&sub, &b, &x) {
                        for (j, &i) in support.iter().enumerate() {
                            if x[j] > RANK_TOL {
                                covered[i] = true;
                            }
                        }
                    }
                }
            }
        }
        if covered.iter().all(|&v| v) {
            return Ok(());
        }
    }
    Err("facet normals do not positively span space at numerical tolerance (unbounded or ill-conditioned input)".into())
}
/// Unpruned enumeration of every support and anchored cyclic order.
pub fn ehz_from_facets(normals: &Vertices, heights: &[f64]) -> Result<ResultData, String> {
    let m = normals.len();
    let d = normals.first().map_or(0, Vec::len);
    if m < d + 1
        || d == 0
        || d % 2 != 0
        || heights.len() != m
        || heights.iter().any(|h| !h.is_finite() || *h <= 0.)
        || normals
            .iter()
            .any(|r| r.len() != d || r.iter().any(|v| !v.is_finite()))
    {
        return Err("invalid facets: require finite even-dimensional normals and positive heights about an interior origin".into());
    }
    // Theorem 1.1 in dimensionless variables t_i = beta_i h_i.
    // Positive row scaling cancels before any rank or KKT calculation.
    let lengths: Vec<f64> = normals
        .iter()
        .map(|r| r.iter().fold(0_f64, |a, &v| a.hypot(v)))
        .collect();
    if lengths.iter().any(|v| !v.is_finite() || *v == 0.) {
        return Err("zero or nonfinite facet normal length".into());
    }
    let unit_heights: Vec<f64> = heights.iter().zip(&lengths).map(|(h, n)| h / n).collect();
    let length_scale = unit_heights.iter().copied().fold(f64::INFINITY, f64::min);
    if length_scale <= 0. || !length_scale.is_finite() {
        return Err("unrepresentable facet length scale".into());
    }
    let p: Vertices = normals
        .iter()
        .zip(&lengths)
        .zip(&unit_heights)
        .map(|((r, n), h)| r.iter().map(|v| (v / n) * (length_scale / h)).collect())
        .collect();
    validate_bounded(&p)?;
    let pairing = Matrix::from_fn(m, m, |i, j| omega(&p[i], &p[j]));
    let mut result = ResultData {
        capacity: 0.,
        minimum: f64::INFINITY,
        order: vec![],
        beta: vec![],
        supports: 0,
        orders: 0,
        singular_kkt: 0,
        status: "floating-point exhaustive theorem-1.1 estimate; not certified",
        maximum: 0.,
        maximizing_order: vec![],
        maximizing_beta: vec![],
        closure_residual: 0.,
        normalization_residual: 0.,
    };
    for size in 2..=m {
        for support in (0..m).combinations(size) {
            result.supports += 1;
            let a = Matrix::from_fn(
                d + 1,
                size,
                |i, j| {
                    if i == d {
                        1.
                    } else {
                        p[support[j]][i]
                    }
                },
            );
            let b = Vector::from_fn(d + 1, |i, _| if i == d { 1. } else { 0. });
            let full_a = a.clone();
            let full_b = b.clone();
            let (a, b, unique) = match prepare_aeq_and_beq(a, b) {
                EqualitySystem::Inconsistent => continue,
                EqualitySystem::Unique(a, b) => {
                    let Some(beta) = solve_nonsingular(&a, &b).and_then(nonnegative) else {
                        continue;
                    };
                    if !feasible(&full_a, &full_b, &beta) {
                        continue;
                    }
                    (a, b, Some(beta))
                }
                EqualitySystem::Multiple(a, b) => (a, b, None),
            };
            // MATLAB perms returns reverse lexicographic order.
            for tail in (1..size).rev().permutations(size - 1) {
                result.orders += 1;
                let route: Vec<usize> = std::iter::once(0).chain(tail).collect();
                let h = Matrix::from_fn(size, size, |i, j| {
                    if i == j {
                        0.
                    } else {
                        let (r, c) = if i > j { (i, j) } else { (j, i) };
                        pairing[(support[route[r]], support[route[c]])]
                    }
                });
                let beta = if let Some(ref beta) = unique {
                    Vector::from_iterator(size, route.iter().map(|&i| beta[i]))
                } else {
                    let objective_scale = h.amax();
                    if objective_scale == 0. {
                        continue;
                    }
                    let rank = a.nrows();
                    let mat = Matrix::from_fn(size + rank, size + rank, |i, j| {
                        match (i < size, j < size) {
                            (true, true) => h[(i, j)] / objective_scale,
                            (true, false) => a[(j - size, route[i])],
                            (false, true) => a[(i - size, route[j])],
                            (false, false) => 0.,
                        }
                    });
                    let rhs = Vector::from_fn(
                        size + rank,
                        |i, _| if i < size { 0. } else { b[i - size] },
                    );
                    let Some(solution) = solve_nonsingular(&mat, &rhs) else {
                        result.singular_kkt += 1;
                        continue;
                    };
                    solution.rows(0, size).into_owned()
                };
                let Some(beta) = nonnegative(beta) else {
                    continue;
                };
                let ordered_a = Matrix::from_fn(d + 1, size, |i, j| full_a[(i, route[j])]);
                if !feasible(&ordered_a, &full_b, &beta) {
                    continue;
                }
                let score = beta.dot(&(&h * &beta)) / 2.;
                if score < result.minimum {
                    result.minimum = score;
                    result.order = route.iter().map(|&i| support[i]).collect();
                    result.beta = route
                        .iter()
                        .enumerate()
                        .map(|(j, &i)| beta[j] / heights[support[i]])
                        .collect();
                }
            }
        }
    }
    if !result.minimum.is_finite() || result.minimum >= 0. {
        return Err("no negative feasible score found at numerical tolerances".into());
    }
    result.minimum /= length_scale * length_scale;
    result.capacity = -0.5 / result.minimum;
    if !result.capacity.is_finite() || result.capacity <= 0. || !result.minimum.is_finite() {
        return Err("capacity/objective outside representable floating-point range".into());
    }
    result.maximum = -result.minimum;
    result.maximizing_order = result.order.iter().rev().copied().collect();
    result.maximizing_beta = result.beta.iter().rev().copied().collect();
    let closure = Vector::<f64>::from_fn(d, |i, _| {
        result
            .order
            .iter()
            .zip(&result.beta)
            .map(|(&j, b)| normals[j][i] * b)
            .sum()
    });
    result.closure_residual = closure.norm();
    result.normalization_residual = (result
        .order
        .iter()
        .zip(&result.beta)
        .map(|(&j, b)| heights[j] * b)
        .sum::<f64>()
        - 1.)
        .abs();
    Ok(result)
}
pub fn ehz_perms(k: &Vertices) -> Result<(Geometry, ResultData), String> {
    let g = geometry(k)?;
    let result = ehz_from_facets(&g.normals, &g.heights)?;
    Ok((g, result))
}
