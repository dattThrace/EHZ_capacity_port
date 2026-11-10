//! Exact checks and certificates for rational H-input.
use super::arithmetic::{bits, matrix, reduce, scalar, vector, Mat, Rat};
use super::{Budget, ExactInput};
use itertools::Itertools;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RowMap {
    pub original: usize,
    #[serde(with = "scalar")]
    pub scale: Rat,
    pub representative: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Removal {
    pub row: usize,
    pub using: Vec<usize>,
    #[serde(with = "vector")]
    pub multipliers: Vec<Rat>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Domain {
    /// Primitive integer rows, one for each retained original row identity.
    #[serde(with = "matrix")]
    pub normals: Mat,
    #[serde(with = "vector")]
    pub heights: Vec<Rat>,
    /// Canonicalization and positive-duplicate identity for every original row.
    pub rows: Vec<RowMap>,
    pub retained: Vec<usize>,
    /// Strictly positive, normalized closure certificate in ORIGINAL coordinates.
    #[serde(with = "vector")]
    pub bounded_weights: Vec<Rat>,
    /// Farkas certificates, in the order inequalities were removed.
    pub removals: Vec<Removal>,
    pub peak_bits: u64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DomainError {
    Invalid(String),
    Incomplete(String),
}

fn unique_tracked(a: &Mat, b: &[Rat], n: usize, peak_bits: &mut u64) -> Option<Vec<Rat>> {
    let reduced = reduce(a, b, n);
    *peak_bits = (*peak_bits).max(reduced.peak_bits);
    if reduced.inconsistent || reduced.pivots.len() != n {
        return None;
    }
    let mut solution = vec![Rat::zero(); n];
    for (row, &pivot) in reduced.pivots.iter().enumerate() {
        solution[pivot] = reduced.rhs[row].clone();
    }
    Some(solution)
}

fn positive_spanning(
    input: &ExactInput,
    budget: &mut Budget,
    peak_bits: &mut u64,
) -> Result<Vec<Rat>, DomainError> {
    let m = input.normals.len();
    let d = input.normals[0].len();
    let span = reduce(&input.normals, &vec![Rat::zero(); m], d);
    *peak_bits = (*peak_bits).max(span.peak_bits);
    if span.pivots.len() != d {
        return Err(DomainError::Invalid(
            "facet normals do not span the ambient space (unbounded input)".into(),
        ));
    }
    // Vertices of {w >= 0: A^T w = 0, sum(w)=1} have at most d+1
    // nonzero entries. Every coordinate which can be positive is positive at
    // some such vertex; summing those vertices yields a strict certificate.
    let mut accumulated = vec![Rat::zero(); m];
    let rhs: Vec<Rat> = (0..=d)
        .map(|i| if i == d { Rat::one() } else { Rat::zero() })
        .collect();
    for size in 2..=(d + 1).min(m) {
        for support in (0..m).combinations(size) {
            budget.tick().map_err(DomainError::Incomplete)?;
            let a: Mat = (0..=d)
                .map(|j| {
                    support
                        .iter()
                        .map(|&i| {
                            if j == d {
                                Rat::one()
                            } else {
                                input.normals[i][j].clone()
                            }
                        })
                        .collect()
                })
                .collect();
            let Some(w) = unique_tracked(&a, &rhs, size, peak_bits) else {
                continue;
            };
            if w.iter().any(Rat::is_negative) {
                continue;
            }
            for (&i, weight) in support.iter().zip(w) {
                accumulated[i] += weight;
                *peak_bits = (*peak_bits).max(bits(&accumulated[i]));
            }
            if accumulated.iter().all(Rat::is_positive) {
                let mut total = Rat::zero();
                for weight in &accumulated {
                    total += weight;
                    *peak_bits = (*peak_bits).max(bits(&total));
                }
                for weight in &mut accumulated {
                    *weight /= &total;
                    *peak_bits = (*peak_bits).max(bits(weight));
                }
                return Ok(accumulated);
            }
        }
    }
    Err(DomainError::Invalid(
        "facet normals do not positively span the ambient space (unbounded input)".into(),
    ))
}

/// A positive scale makes the full inequality a primitive integer row.
fn canonical_row(normal: &[Rat], height: &Rat, peak_bits: &mut u64) -> (Vec<Rat>, Rat, Rat) {
    // Observe each value constructed by this module, including denominator
    // clearing before cancellation. BigInt gcd internals and num-rational's
    // internal normalization temporaries are outside this measurement.
    let mut denominator = BigInt::one();
    for value in normal.iter().chain(std::iter::once(height)) {
        let gcd = denominator.gcd(value.denom());
        *peak_bits = (*peak_bits).max(gcd.bits());
        let quotient = &denominator / gcd;
        *peak_bits = (*peak_bits).max(quotient.bits());
        denominator = quotient * value.denom();
        *peak_bits = (*peak_bits).max(denominator.bits());
    }
    let mut divisor = BigInt::zero();
    for value in normal.iter().chain(std::iter::once(height)) {
        let factor = &denominator / value.denom();
        *peak_bits = (*peak_bits).max(factor.bits());
        let numerator = (value.numer() * factor).abs();
        *peak_bits = (*peak_bits).max(numerator.bits());
        divisor = divisor.gcd(&numerator);
        *peak_bits = (*peak_bits).max(divisor.bits());
    }
    let scale = Rat::new(denominator, divisor);
    *peak_bits = (*peak_bits).max(bits(&scale));
    let normal: Vec<_> = normal
        .iter()
        .map(|value| {
            let scaled = value * &scale;
            *peak_bits = (*peak_bits).max(bits(&scaled));
            scaled
        })
        .collect();
    let height = height * &scale;
    *peak_bits = (*peak_bits).max(bits(&height));
    (normal, height, scale)
}

fn redundant(
    input: &ExactInput,
    row: usize,
    current: &[usize],
    budget: &mut Budget,
    peak_bits: &mut u64,
) -> Result<Option<Removal>, DomainError> {
    let d = input.normals[0].len();
    let others: Vec<_> = current.iter().copied().filter(|&i| i != row).collect();
    // Minimize h^T lambda over A^T lambda=a_row, lambda>=0. Positive
    // heights make sublevel sets compact. A feasible optimum has an
    // independent support of size <=d, so this search is complete for the
    // Farkas certificate h^T lambda<=h_row without an LP dependency.
    for size in 1..=d.min(others.len()) {
        for using in others.iter().copied().combinations(size) {
            budget.tick().map_err(DomainError::Incomplete)?;
            let a: Mat = (0..d)
                .map(|j| using.iter().map(|&i| input.normals[i][j].clone()).collect())
                .collect();
            let Some(multipliers) = unique_tracked(&a, &input.normals[row], size, peak_bits) else {
                continue;
            };
            if multipliers.iter().any(Rat::is_negative) {
                continue;
            }
            let mut height = Rat::zero();
            for (&i, lambda) in using.iter().zip(&multipliers) {
                let product = &input.heights[i] * lambda;
                *peak_bits = (*peak_bits).max(bits(&product));
                height += product;
                *peak_bits = (*peak_bits).max(bits(&height));
            }
            if height <= input.heights[row] {
                return Ok(Some(Removal {
                    row,
                    using,
                    multipliers,
                }));
            }
        }
    }
    Ok(None)
}

pub fn prepare(input: &ExactInput, budget: &mut Budget) -> Result<Domain, DomainError> {
    let m = input.normals.len();
    let d = input.normals.first().map_or(0, Vec::len);
    if d == 0
        || d % 2 != 0
        || m < d + 1
        || input.heights.len() != m
        || input.heights.iter().any(|h| h <= &Rat::zero())
        || input
            .normals
            .iter()
            .any(|row| row.len() != d || row.iter().all(Rat::is_zero))
    {
        return Err(DomainError::Invalid(
            "invalid exact facets: require nonzero equal even-dimensional normals and positive heights about an interior origin".into(),
        ));
    }
    let mut peak_bits = input
        .normals
        .iter()
        .flatten()
        .chain(&input.heights)
        .map(bits)
        .max()
        .unwrap_or(0);
    let bounded_weights = positive_spanning(input, budget, &mut peak_bits)?;
    let mut canonical_normals: Mat = Vec::with_capacity(m);
    let mut canonical_heights = Vec::with_capacity(m);
    let mut rows = Vec::with_capacity(m);
    let mut retained = vec![];
    for i in 0..m {
        let (normal, height, scale) =
            canonical_row(&input.normals[i], &input.heights[i], &mut peak_bits);
        let representative = (0..i)
            .find(|&j| canonical_normals[j] == normal && canonical_heights[j] == height)
            .unwrap_or(i);
        if representative == i {
            retained.push(i);
        }
        canonical_normals.push(normal);
        canonical_heights.push(height);
        rows.push(RowMap {
            original: i,
            scale,
            representative,
        });
    }
    let mut removals = vec![];
    // Sequential removal is essential: each proof only uses inequalities
    // still present when the certified deletion is applied.
    for row in retained.clone() {
        if let Some(removal) = redundant(input, row, &retained, budget, &mut peak_bits)? {
            retained.retain(|&i| i != row);
            removals.push(removal);
        }
    }
    let normals = retained
        .iter()
        .map(|&i| canonical_normals[i].clone())
        .collect();
    let heights = retained
        .iter()
        .map(|&i| canonical_heights[i].clone())
        .collect();
    Ok(Domain {
        normals,
        heights,
        rows,
        retained,
        bounded_weights,
        removals,
        peak_bits,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::exact::arithmetic::parse_rat;
    use crate::exact::SearchOptions;

    fn input(normals: &[&[&str]], heights: &[&str]) -> ExactInput {
        ExactInput {
            schema: "ehz-exact-h-v1".into(),
            normals: normals
                .iter()
                .map(|row| row.iter().map(|v| parse_rat(v).unwrap()).collect())
                .collect(),
            heights: heights.iter().map(|v| parse_rat(v).unwrap()).collect(),
        }
    }

    fn square() -> ExactInput {
        input(
            &[&["1", "0"], &["-1", "0"], &["0", "1"], &["0", "-1"]],
            &["1", "1", "1", "1"],
        )
    }

    fn check_certificates(input: &ExactInput, domain: &Domain) {
        assert!(domain.bounded_weights.iter().all(Rat::is_positive));
        assert_eq!(
            domain.bounded_weights.iter().cloned().sum::<Rat>(),
            Rat::one()
        );
        for j in 0..input.normals[0].len() {
            assert!(input
                .normals
                .iter()
                .zip(&domain.bounded_weights)
                .map(|(a, w)| &a[j] * w)
                .sum::<Rat>()
                .is_zero());
        }
        let mut current: Vec<usize> = domain
            .rows
            .iter()
            .filter(|row| row.original == row.representative)
            .map(|row| row.original)
            .collect();
        for removal in &domain.removals {
            assert!(current.contains(&removal.row));
            assert!(removal.using.iter().all(|i| current.contains(i)));
            assert!(!removal.using.contains(&removal.row));
            assert!(removal.multipliers.iter().all(|v| v >= &Rat::zero()));
            for j in 0..input.normals[0].len() {
                let lhs: Rat = removal
                    .using
                    .iter()
                    .zip(&removal.multipliers)
                    .map(|(&i, lambda)| &input.normals[i][j] * lambda)
                    .sum();
                assert_eq!(lhs, input.normals[removal.row][j]);
            }
            let height: Rat = removal
                .using
                .iter()
                .zip(&removal.multipliers)
                .map(|(&i, lambda)| &input.heights[i] * lambda)
                .sum();
            assert!(height <= input.heights[removal.row]);
            current.retain(|&i| i != removal.row);
        }
        assert_eq!(current, domain.retained);
    }

    #[test]
    fn boundedness_certificate_covers_original_duplicate_rows() {
        let mut input = square();
        input
            .normals
            .push(vec![parse_rat("2/3").unwrap(), Rat::zero()]);
        input.heights.push(parse_rat("2/3").unwrap());
        let domain = prepare(&input, &mut Budget::new(&SearchOptions::default())).unwrap();
        assert_eq!(domain.retained, vec![0, 1, 2, 3]);
        assert_eq!(domain.rows[4].representative, 0);
        assert_eq!(domain.rows[4].scale, parse_rat("3/2").unwrap());
        assert!(domain.removals.is_empty());
        check_certificates(&input, &domain);
    }

    #[test]
    fn certified_redundancy_includes_tight_vertex_and_weak_inequalities() {
        let mut input = square();
        for (normal, height) in [(["1", "0"], "3"), (["1", "1"], "2")] {
            input
                .normals
                .push(normal.iter().map(|v| parse_rat(v).unwrap()).collect());
            input.heights.push(parse_rat(height).unwrap());
        }
        let domain = prepare(&input, &mut Budget::new(&SearchOptions::default())).unwrap();
        assert_eq!(domain.retained, vec![0, 1, 2, 3]);
        assert_eq!(domain.removals.len(), 2);
        check_certificates(&input, &domain);
    }

    #[test]
    fn genuine_truncation_is_retained() {
        let mut input = square();
        input.normals.push(vec![Rat::one(), Rat::one()]);
        input.heights.push(Rat::one());
        let domain = prepare(&input, &mut Budget::new(&SearchOptions::default())).unwrap();
        assert_eq!(domain.retained, vec![0, 1, 2, 3, 4]);
        check_certificates(&input, &domain);
    }

    #[test]
    fn unbounded_and_malformed_domains_are_invalid() {
        let unbounded = input(&[&["1", "0"], &["-1", "0"], &["0", "1"]], &["1", "1", "1"]);
        assert!(matches!(
            prepare(&unbounded, &mut Budget::new(&SearchOptions::default())),
            Err(DomainError::Invalid(_))
        ));
        let mut malformed = square();
        malformed.heights[0] = Rat::zero();
        assert!(matches!(
            prepare(&malformed, &mut Budget::new(&SearchOptions::default())),
            Err(DomainError::Invalid(_))
        ));
        let mut zero_normal = square();
        zero_normal.normals[0].fill(Rat::zero());
        assert!(matches!(
            prepare(&zero_normal, &mut Budget::new(&SearchOptions::default())),
            Err(DomainError::Invalid(_))
        ));
    }

    #[test]
    fn exhausted_domain_budget_is_incomplete() {
        let options = SearchOptions {
            max_steps: Some(0),
            ..SearchOptions::default()
        };
        assert!(matches!(
            prepare(&square(), &mut Budget::new(&options)),
            Err(DomainError::Incomplete(_))
        ));
    }
}
