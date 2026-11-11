//! Exact reference solver for small rational H-polytopes.
//!
//! Every support and anchored cyclic order is covered. A singular stationary
//! system is omitted only because its optimum is attained on a smaller support.
//! No floating-point arithmetic, rank tolerance, or geometric adapter is used.
pub mod arithmetic;
pub mod domain;

use arithmetic::{bits, matrix, reduce, scalar, vector, Mat, Rat};
use domain::{Domain, DomainError};
use itertools::Itertools;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::time::Instant;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ExactInput {
    pub schema: String,
    #[serde(with = "matrix")]
    pub normals: Mat,
    #[serde(with = "vector")]
    pub heights: Vec<Rat>,
}

#[derive(Clone, Debug, Default)]
pub struct SearchOptions {
    /// Work items across domain checks, supports, and orders; None is unlimited.
    pub max_steps: Option<u64>,
    /// Cooperative deadline, checked between work items, not inside elimination.
    pub timeout_ms: Option<u64>,
}

pub struct Budget {
    started: Instant,
    max_steps: Option<u64>,
    timeout_ms: Option<u64>,
    pub steps: u64,
}

impl Budget {
    pub fn new(options: &SearchOptions) -> Self {
        Self {
            started: Instant::now(),
            max_steps: options.max_steps,
            timeout_ms: options.timeout_ms,
            steps: 0,
        }
    }

    pub fn tick(&mut self) -> Result<(), String> {
        if self.max_steps.is_some_and(|limit| self.steps >= limit) {
            return Err("work budget exhausted".into());
        }
        if self
            .timeout_ms
            .is_some_and(|limit| self.started.elapsed().as_millis() >= u128::from(limit))
        {
            return Err("cooperative timeout reached".into());
        }
        self.steps = self.steps.saturating_add(1);
        Ok(())
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Witness {
    /// Original input row identities, in ordered objective order.
    pub order: Vec<usize>,
    /// Weights in ORIGINAL, unnormalized input coordinates, parallel to order.
    #[serde(with = "vector")]
    pub weights: Vec<Rat>,
    #[serde(with = "scalar")]
    pub objective: Rat,
    /// A positive feasible objective alone certifies only this upper bound.
    #[serde(with = "scalar")]
    pub capacity_upper_bound: Rat,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Completion {
    pub complete: bool,
    pub reason: String,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct CapacityInterval {
    #[serde(with = "scalar")]
    pub lower: Rat,
    #[serde(with = "scalar")]
    pub upper: Rat,
    #[serde(with = "scalar")]
    pub objective_upper_bound: Rat,
    pub method: String,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct SearchStats {
    pub dimension: usize,
    pub original_facets: usize,
    pub retained_facets: usize,
    pub domain_steps: u64,
    pub steps: u64,
    pub supports: u64,
    pub orders: u64,
    pub inconsistent_equalities: u64,
    pub unique_equalities: u64,
    pub free_equalities: u64,
    pub inconsistent_stationary: u64,
    pub singular_stationary: u64,
    pub feasible_candidates: u64,
    pub input_bits: u64,
    /// Observed rational components, including elimination products; excludes
    /// transient BigInt internals inside rational normalization.
    pub peak_rational_bits: u64,
    pub elapsed_ms: u128,
    /// Linux process high-water RSS; includes runtime and allocator overhead.
    pub peak_rss_kib: Option<u64>,
}

pub mod optional_scalar {
    use super::*;
    pub fn serialize<S: serde::Serializer>(value: &Option<Rat>, s: S) -> Result<S::Ok, S::Error> {
        value.as_ref().map(arithmetic::canonical).serialize(s)
    }
    pub fn deserialize<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Option<Rat>, D::Error> {
        Option::<String>::deserialize(d)?
            .map(|value| arithmetic::parse_rat(&value).map_err(serde::de::Error::custom))
            .transpose()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExactResult {
    pub schema: String,
    pub status: String,
    /// Canonical semantic input identity, independent of source JSON formatting.
    pub input: ExactInput,
    pub domain: Option<Domain>,
    pub witness: Option<Witness>,
    #[serde(with = "optional_scalar")]
    pub capacity: Option<Rat>,
    pub interval: Option<CapacityInterval>,
    pub completion: Completion,
    pub stats: SearchStats,
    pub error: Option<String>,
}

pub fn omega(a: &[Rat], b: &[Rat]) -> Rat {
    omega_observed(a, b, &mut 0)
}

fn observed(value: Rat, peak: &mut u64) -> Rat {
    *peak = (*peak).max(bits(&value));
    value
}

fn omega_observed(a: &[Rat], b: &[Rat], peak: &mut u64) -> Rat {
    let n = a.len() / 2;
    let mut total = Rat::zero();
    for i in 0..n {
        let left = observed(&a[i] * &b[n + i], peak);
        let right = observed(&a[n + i] * &b[i], peak);
        let difference = observed(left - right, peak);
        total = observed(total + difference, peak);
    }
    total
}

/// Recheck a candidate directly against the ORIGINAL input and rescore the
/// ordered sum. This uses no elimination or stationary-system data. It is an
/// arithmetic residual check; establishing the domain and global optimality
/// requires the separate certificate verifier.
pub fn check_witness(input: &ExactInput, witness: &Witness) -> Result<(), String> {
    check_witness_observed(input, witness, &mut 0)
}

fn check_witness_observed(
    input: &ExactInput,
    witness: &Witness,
    peak: &mut u64,
) -> Result<(), String> {
    let d = input.normals.first().map_or(0, Vec::len);
    if d == 0
        || input.heights.len() != input.normals.len()
        || input.normals.iter().any(|row| row.len() != d)
        || witness.order.len() != witness.weights.len()
        || witness.order.is_empty()
        || witness.order.iter().any(|&i| i >= input.normals.len())
        || witness.order.iter().unique().count() != witness.order.len()
        || witness.weights.iter().any(Rat::is_negative)
    {
        return Err("malformed exact witness".into());
    }
    for j in 0..d {
        let mut residual = Rat::zero();
        for (&i, weight) in witness.order.iter().zip(&witness.weights) {
            let term = observed(&input.normals[i][j] * weight, peak);
            residual = observed(residual + term, peak);
        }
        if !residual.is_zero() {
            return Err("exact witness fails original closure".into());
        }
    }
    let mut mass = Rat::zero();
    for (&i, weight) in witness.order.iter().zip(&witness.weights) {
        let term = observed(&input.heights[i] * weight, peak);
        mass = observed(mass + term, peak);
    }
    if mass != Rat::one() {
        return Err("exact witness fails original normalization".into());
    }
    let score = ordered_score(input, &witness.order, &witness.weights, peak);
    let product = observed(&score * &witness.capacity_upper_bound, peak);
    let reciprocal_check = observed(product * Rat::from_integer(2.into()), peak);
    if !score.is_positive() || score != witness.objective || reciprocal_check != Rat::one() {
        return Err("exact witness fails objective or reciprocal check".into());
    }
    Ok(())
}

fn ordered_score(input: &ExactInput, order: &[usize], weights: &[Rat], peak: &mut u64) -> Rat {
    let mut score = Rat::zero();
    for i in 0..order.len() {
        for j in 0..i {
            let pairing = omega_observed(&input.normals[order[i]], &input.normals[order[j]], peak);
            let product = observed(&weights[i] * &weights[j], peak);
            let term = observed(product * pairing, peak);
            score = observed(score + term, peak);
        }
    }
    score
}

/// Optionally turn a partial positive witness into a certified interval using
/// a conservative envelope. With t_i=beta_i h_i and M=max |omega(a_i/h_i,a_j/h_j)|,
/// Q <= M sum_{i>j} t_i t_j <= M/2. No search region is discarded by this bound.
pub fn certify_interval(result: &mut ExactResult) -> Result<(), String> {
    if result.status != "exact-witness/incomplete" {
        return Ok(());
    }
    let witness = result.witness.as_ref().ok_or("missing positive witness")?;
    let peak = &mut result.stats.peak_rational_bits;
    check_witness_observed(&result.input, witness, peak)?;
    if result.domain.is_none() || result.completion.complete || result.capacity.is_some() {
        return Err("interval requires a certified domain and incomplete search".into());
    }
    let p: Mat = result
        .input
        .normals
        .iter()
        .zip(&result.input.heights)
        .map(|(row, h)| row.iter().map(|v| observed(v / h, peak)).collect())
        .collect();
    let mut maximum = Rat::zero();
    for i in 0..p.len() {
        for j in 0..i {
            maximum = maximum.max(omega_observed(&p[i], &p[j], peak).abs());
        }
    }
    let upper = observed(&maximum / Rat::from_integer(2.into()), peak);
    if upper < witness.objective || !upper.is_positive() {
        return Err("pairing envelope contradicts the positive witness".into());
    }
    result.interval = Some(CapacityInterval {
        lower: observed(Rat::one() / maximum, peak),
        upper: witness.capacity_upper_bound.clone(),
        objective_upper_bound: upper,
        method: "pairing-envelope-v1".into(),
    });
    result.status = "certified-interval".into();
    Ok(())
}

fn resident_peak() -> Option<u64> {
    std::fs::read_to_string("/proc/self/status")
        .ok()?
        .lines()
        .find(|line| line.starts_with("VmHWM:"))?
        .split_whitespace()
        .nth(1)?
        .parse()
        .ok()
}

/// Reference exhaustive exact search. Complete status is emitted only after
/// domain certification and the entire support/order traversal finishes.
pub fn solve(input: &ExactInput, options: &SearchOptions) -> ExactResult {
    let mut budget = Budget::new(options);
    let input_bits = input
        .normals
        .iter()
        .flatten()
        .chain(&input.heights)
        .map(bits)
        .max()
        .unwrap_or(0);
    let mut result = ExactResult {
        schema: "ehz-exact-result-v1".into(),
        status: "incomplete".into(),
        input: input.clone(),
        domain: None,
        witness: None,
        capacity: None,
        interval: None,
        completion: Completion {
            complete: false,
            reason: "search not completed".into(),
        },
        stats: SearchStats {
            dimension: input.normals.first().map_or(0, Vec::len),
            original_facets: input.normals.len(),
            input_bits,
            peak_rational_bits: input_bits,
            ..SearchStats::default()
        },
        error: None,
    };
    if input.schema != "ehz-exact-h-v1" {
        result.status = "invalid-input".into();
        result.error = Some("unsupported exact H-input schema".into());
    } else {
        match domain::prepare(input, &mut budget) {
            Err(DomainError::Invalid(error)) => {
                result.status = "invalid-input".into();
                result.error = Some(error);
            }
            Err(DomainError::Incomplete(reason)) => result.completion.reason = reason,
            Ok(domain) => {
                result.stats.domain_steps = budget.steps;
                result.stats.retained_facets = domain.retained.len();
                result.stats.peak_rational_bits =
                    result.stats.peak_rational_bits.max(domain.peak_bits);
                let outcome = enumerate(input, &domain, &mut budget, &mut result);
                result.domain = Some(domain);
                match outcome {
                    Ok(()) if result.witness.is_some() => {
                        result.status = "exact-complete".into();
                        result.capacity = result
                            .witness
                            .as_ref()
                            .map(|w| w.capacity_upper_bound.clone());
                        result.completion = Completion {
                            complete: true,
                            reason: "all supports and anchored cyclic orders exhausted; smaller supports cover exactly singular stationary extrema".into(),
                        };
                    }
                    Ok(()) => {
                        result.error = Some("internal consistency failure: no positive objective on a certified bounded full-dimensional body".into());
                        result.completion.reason = "no verified positive witness".into();
                    }
                    Err(reason) => {
                        result.completion.reason = reason;
                        if result.witness.is_some() {
                            result.status = "exact-witness/incomplete".into();
                        }
                    }
                }
            }
        }
    }
    if result.status == "invalid-input" {
        result.completion.reason = "input rejected before capacity search".into();
    }
    result.stats.steps = budget.steps;
    if result.domain.is_none() {
        result.stats.domain_steps = budget.steps;
    }
    result.stats.elapsed_ms = budget.started.elapsed().as_millis();
    result.stats.peak_rss_kib = resident_peak();
    result
}

fn enumerate(
    input: &ExactInput,
    domain: &Domain,
    budget: &mut Budget,
    result: &mut ExactResult,
) -> Result<(), String> {
    let m = domain.normals.len();
    let d = domain.normals[0].len();
    let pairings: Mat = (0..m)
        .map(|i| {
            (0..m)
                .map(|j| {
                    omega_observed(
                        &domain.normals[i],
                        &domain.normals[j],
                        &mut result.stats.peak_rational_bits,
                    )
                })
                .collect()
        })
        .collect();
    for value in pairings.iter().flatten() {
        result.stats.peak_rational_bits = result.stats.peak_rational_bits.max(bits(value));
    }
    for size in 2..=m {
        for support in (0..m).combinations(size) {
            budget.tick()?;
            result.stats.supports += 1;
            let a: Mat = (0..=d)
                .map(|j| {
                    support
                        .iter()
                        .map(|&i| {
                            if j == d {
                                domain.heights[i].clone()
                            } else {
                                domain.normals[i][j].clone()
                            }
                        })
                        .collect()
                })
                .collect();
            let rhs: Vec<Rat> = (0..=d)
                .map(|j| if j == d { Rat::one() } else { Rat::zero() })
                .collect();
            // Cache this exact reduction for every order on the same support.
            let equality = reduce(&a, &rhs, size);
            result.stats.peak_rational_bits =
                result.stats.peak_rational_bits.max(equality.peak_bits);
            if equality.inconsistent {
                result.stats.inconsistent_equalities += 1;
                continue;
            }
            let unique = equality.pivots.len() == size;
            if unique {
                result.stats.unique_equalities += 1;
                if equality.rhs.iter().any(Rat::is_negative) {
                    continue;
                }
            } else {
                result.stats.free_equalities += 1;
            }
            for tail in (1..size).permutations(size - 1) {
                budget.tick()?;
                result.stats.orders += 1;
                let route: Vec<usize> = std::iter::once(0).chain(tail).collect();
                let mut h = vec![vec![Rat::zero(); size]; size];
                for i in 0..size {
                    for j in 0..i {
                        let (u, v) = (route[i], route[j]);
                        h[u][v] = pairings[support[u]][support[v]].clone();
                        h[v][u] = h[u][v].clone();
                    }
                }
                let beta = if unique {
                    equality.rhs.clone()
                } else {
                    let rank = equality.pivots.len();
                    let kkt: Mat = (0..size + rank)
                        .map(|i| {
                            (0..size + rank)
                                .map(|j| match (i < size, j < size) {
                                    (true, true) => h[i][j].clone(),
                                    (true, false) => equality.rows[j - size][i].clone(),
                                    (false, true) => equality.rows[i - size][j].clone(),
                                    (false, false) => Rat::zero(),
                                })
                                .collect()
                        })
                        .collect();
                    let kkt_rhs: Vec<Rat> = (0..size + rank)
                        .map(|i| {
                            if i < size {
                                Rat::zero()
                            } else {
                                equality.rhs[i - size].clone()
                            }
                        })
                        .collect();
                    let stationary = reduce(&kkt, &kkt_rhs, size + rank);
                    result.stats.peak_rational_bits =
                        result.stats.peak_rational_bits.max(stationary.peak_bits);
                    if stationary.inconsistent {
                        result.stats.inconsistent_stationary += 1;
                        continue;
                    }
                    if stationary.pivots.len() != size + rank {
                        result.stats.singular_stationary += 1;
                        continue;
                    }
                    stationary.rhs[..size].to_vec()
                };
                if beta.iter().any(Rat::is_negative) {
                    continue;
                }
                // Check every original equality of this support, including all
                // rows discarded by exact row reduction.
                if a.iter().zip(&rhs).any(|(row, target)| {
                    row.iter().zip(&beta).map(|(x, y)| x * y).sum::<Rat>() != *target
                }) {
                    return Err("internal exact equality residual failure".into());
                }
                let order: Vec<usize> =
                    route.iter().map(|&i| domain.retained[support[i]]).collect();
                let weights: Vec<Rat> = route
                    .iter()
                    .map(|&i| {
                        observed(
                            &beta[i] * &domain.rows[domain.retained[support[i]]].scale,
                            &mut result.stats.peak_rational_bits,
                        )
                    })
                    .collect();
                let objective = ordered_score(
                    input,
                    &order,
                    &weights,
                    &mut result.stats.peak_rational_bits,
                );
                // Independently compare original ordered sum with the KKT
                // quadratic in normalized support coordinates.
                let mut quadratic = Rat::zero();
                for i in 0..size {
                    for j in 0..size {
                        let term =
                            observed(&beta[i] * &h[i][j], &mut result.stats.peak_rational_bits);
                        let term = observed(term * &beta[j], &mut result.stats.peak_rational_bits);
                        quadratic =
                            observed(quadratic + term, &mut result.stats.peak_rational_bits);
                    }
                }
                if quadratic != &objective * Rat::from_integer(2.into()) {
                    return Err("internal exact objective mapping failure".into());
                }
                result.stats.feasible_candidates += 1;
                result.stats.peak_rational_bits = result
                    .stats
                    .peak_rational_bits
                    .max(bits(&objective))
                    .max(bits(&quadratic));
                if objective.is_positive()
                    && result
                        .witness
                        .as_ref()
                        .is_none_or(|best| objective > best.objective)
                {
                    let witness = Witness {
                        order,
                        weights,
                        capacity_upper_bound: observed(
                            Rat::one() / (Rat::from_integer(2.into()) * &objective),
                            &mut result.stats.peak_rational_bits,
                        ),
                        objective,
                    };
                    check_witness_observed(input, &witness, &mut result.stats.peak_rational_bits)?;
                    result.witness = Some(witness);
                }
            }
        }
    }
    Ok(())
}
