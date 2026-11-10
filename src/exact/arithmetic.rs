//! Canonical rational input and transparent exact row reduction.
//!
//! JSON rationals are strings or pairs of integer strings. JSON numbers are
//! deliberately excluded: accepting them through a floating decoder would lose
//! the exact input before the solver sees it.
use num_bigint::BigInt;
use num_traits::{One, Zero};
use serde::{Deserialize, Deserializer, Serializer};
use std::str::FromStr;

pub type Rat = num_rational::BigRational;
pub type Mat = Vec<Vec<Rat>>;

fn integer(text: &str) -> Result<BigInt, String> {
    let digits = text.strip_prefix(['+', '-']).unwrap_or(text);
    if digits.is_empty() || !digits.bytes().all(|c| c.is_ascii_digit()) {
        return Err(format!("invalid integer: {text:?}"));
    }
    BigInt::from_str(text).map_err(|_| format!("invalid integer: {text:?}"))
}

/// Parse an integer, exact decimal (optionally scientific), or integer fraction.
/// Decimal exponents are limited to one million to bound accidental expansion;
/// ordinary integer and fraction inputs have no fixed precision limit.
pub fn parse_rat(text: &str) -> Result<Rat, String> {
    if text.contains('/') {
        let (numerator, denominator) = text
            .split_once('/')
            .ok_or_else(|| "invalid fraction".to_string())?;
        let numerator = integer(numerator)?;
        let denominator = integer(denominator)?;
        if denominator.is_zero() {
            return Err("zero denominator".into());
        }
        return Ok(Rat::new(numerator, denominator));
    }
    let mut parts = text.split(['e', 'E']);
    let mantissa = parts.next().unwrap_or("");
    let exponent = match parts.next() {
        None => 0_i64,
        Some(value) => {
            integer(value)?;
            value
                .parse::<i64>()
                .map_err(|_| "decimal exponent too large".to_string())?
        }
    };
    if parts.next().is_some() {
        return Err("invalid decimal exponent".into());
    }
    let (negative, unsigned) = match mantissa.as_bytes().first() {
        Some(b'-') => (true, &mantissa[1..]),
        Some(b'+') => (false, &mantissa[1..]),
        _ => (false, mantissa),
    };
    let mut decimal_parts = unsigned.split('.');
    let integral = decimal_parts.next().unwrap_or("");
    let fractional = decimal_parts.next().unwrap_or("");
    if decimal_parts.next().is_some()
        || (integral.is_empty() && fractional.is_empty())
        || !integral.bytes().all(|c| c.is_ascii_digit())
        || !fractional.bytes().all(|c| c.is_ascii_digit())
    {
        return Err(format!("invalid rational: {text:?}"));
    }
    let places = i64::try_from(fractional.len()).map_err(|_| "decimal too long")?;
    let shift = exponent
        .checked_sub(places)
        .ok_or("decimal exponent too large")?;
    if shift.unsigned_abs() > 1_000_000 {
        return Err("decimal exponent too large".into());
    }
    let digits = format!("{integral}{fractional}");
    let mut numerator = BigInt::from_str(&digits).map_err(|_| "invalid decimal")?;
    if negative {
        numerator = -numerator;
    }
    let power = BigInt::from(10_u8).pow(shift.unsigned_abs() as u32);
    Ok(if shift >= 0 {
        Rat::from_integer(numerator * power)
    } else {
        Rat::new(numerator, power)
    })
}

pub fn canonical(value: &Rat) -> String {
    if value.denom().is_one() {
        value.numer().to_string()
    } else {
        format!("{}/{}", value.numer(), value.denom())
    }
}

/// Largest magnitude bit length of a rational's numerator and denominator.
pub fn bits(value: &Rat) -> u64 {
    value.numer().bits().max(value.denom().bits())
}

#[derive(Deserialize)]
#[serde(untagged)]
enum Encoded {
    Text(String),
    Pair([String; 2]),
}

impl Encoded {
    fn rational(self) -> Result<Rat, String> {
        match self {
            Self::Text(text) => parse_rat(&text),
            Self::Pair([numerator, denominator]) => {
                let numerator = integer(&numerator)?;
                let denominator = integer(&denominator)?;
                if denominator.is_zero() {
                    Err("zero denominator".into())
                } else {
                    Ok(Rat::new(numerator, denominator))
                }
            }
        }
    }
}

pub mod scalar {
    use super::*;

    pub fn serialize<S: Serializer>(value: &Rat, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&canonical(value))
    }

    pub fn deserialize<'de, D: Deserializer<'de>>(deserializer: D) -> Result<Rat, D::Error> {
        Encoded::deserialize(deserializer)?
            .rational()
            .map_err(serde::de::Error::custom)
    }
}

pub mod vector {
    use super::*;
    use serde::Serialize;

    pub fn serialize<S: Serializer>(values: &[Rat], serializer: S) -> Result<S::Ok, S::Error> {
        values
            .iter()
            .map(canonical)
            .collect::<Vec<_>>()
            .serialize(serializer)
    }

    pub fn deserialize<'de, D: Deserializer<'de>>(deserializer: D) -> Result<Vec<Rat>, D::Error> {
        Vec::<Encoded>::deserialize(deserializer)?
            .into_iter()
            .map(|value| value.rational().map_err(serde::de::Error::custom))
            .collect()
    }
}

pub mod matrix {
    use super::*;
    use serde::Serialize;

    pub fn serialize<S: Serializer>(values: &Mat, serializer: S) -> Result<S::Ok, S::Error> {
        values
            .iter()
            .map(|row| row.iter().map(canonical).collect::<Vec<_>>())
            .collect::<Vec<_>>()
            .serialize(serializer)
    }

    pub fn deserialize<'de, D: Deserializer<'de>>(deserializer: D) -> Result<Mat, D::Error> {
        Vec::<Vec<Encoded>>::deserialize(deserializer)?
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|value| value.rational().map_err(serde::de::Error::custom))
                    .collect()
            })
            .collect()
    }
}

#[derive(Debug, Clone)]
pub struct Reduction {
    /// Independent reduced rows. Pivot columns are in increasing order.
    pub rows: Mat,
    pub rhs: Vec<Rat>,
    pub pivots: Vec<usize>,
    pub inconsistent: bool,
    /// Includes products before cancellation, as well as stored row entries.
    pub peak_bits: u64,
}

/// Exact Gauss-Jordan elimination with row pivoting and explicit free columns.
/// Shape errors are represented as inconsistent systems, never truncated rows.
pub fn reduce(a: &Mat, b: &[Rat], variables: usize) -> Reduction {
    let mut result = Reduction {
        rows: a.clone(),
        rhs: b.to_vec(),
        pivots: vec![],
        inconsistent: false,
        peak_bits: a.iter().flatten().chain(b).map(bits).max().unwrap_or(0),
    };
    if a.len() != b.len() || a.iter().any(|row| row.len() != variables) {
        result.inconsistent = true;
        return result;
    }
    for column in 0..variables {
        let target = result.pivots.len();
        let Some(pivot) = (target..result.rows.len()).find(|&i| !result.rows[i][column].is_zero())
        else {
            continue;
        };
        result.rows.swap(target, pivot);
        result.rhs.swap(target, pivot);
        let scale = result.rows[target][column].clone();
        for value in &mut result.rows[target] {
            *value /= &scale;
            result.peak_bits = result.peak_bits.max(bits(value));
        }
        result.rhs[target] /= &scale;
        result.peak_bits = result.peak_bits.max(bits(&result.rhs[target]));
        for row in 0..result.rows.len() {
            if row == target || result.rows[row][column].is_zero() {
                continue;
            }
            let multiplier = result.rows[row][column].clone();
            for col in 0..variables {
                let product = &multiplier * &result.rows[target][col];
                result.peak_bits = result.peak_bits.max(bits(&product));
                result.rows[row][col] -= product;
                result.peak_bits = result.peak_bits.max(bits(&result.rows[row][col]));
            }
            let product = multiplier * &result.rhs[target];
            result.peak_bits = result.peak_bits.max(bits(&product));
            result.rhs[row] -= product;
            result.peak_bits = result.peak_bits.max(bits(&result.rhs[row]));
        }
        result.pivots.push(column);
    }
    let rank = result.pivots.len();
    result.inconsistent = result.rhs.iter().skip(rank).any(|rhs| !rhs.is_zero());
    result.rows.truncate(rank);
    result.rhs.truncate(rank);
    result
}

pub fn unique(a: &Mat, b: &[Rat], variables: usize) -> Option<Vec<Rat>> {
    let reduced = reduce(a, b, variables);
    if reduced.inconsistent || reduced.pivots.len() != variables {
        None
    } else {
        Some(reduced.rhs)
    }
}

pub fn rank(a: &Mat, variables: usize) -> usize {
    reduce(a, &vec![Rat::zero(); a.len()], variables)
        .pivots
        .len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde::Serialize;

    fn r(value: &str) -> Rat {
        parse_rat(value).unwrap()
    }

    fn m(rows: &[&[&str]]) -> Mat {
        rows.iter()
            .map(|row| row.iter().map(|s| r(s)).collect())
            .collect()
    }

    #[test]
    fn exact_syntax_and_canonicalization() {
        for (input, expected) in [
            ("-24/-36", "2/3"),
            ("+2.500", "5/2"),
            (".125", "1/8"),
            ("5.", "5"),
            ("12.5e-3", "1/80"),
            ("-1E+3", "-1000"),
            ("0/-7", "0"),
            ("9007199254740993", "9007199254740993"),
            (
                "123456789012345678901234567890123456789012345678901234567891/2",
                "123456789012345678901234567890123456789012345678901234567891/2",
            ),
        ] {
            assert_eq!(canonical(&r(input)), expected);
            assert_eq!(r(input), r(&canonical(&r(input))));
        }
        for input in [
            "",
            ".",
            "+",
            "1/0",
            "1/2/3",
            "1.5/2",
            " 2",
            "1 ",
            "NaN",
            "1e",
            "1ee2",
            "1..2",
            "--1",
            "1e1000001",
        ] {
            assert!(parse_rat(input).is_err(), "accepted {input:?}");
        }
    }

    #[test]
    fn serde_never_accepts_machine_numbers() {
        #[derive(Debug, Serialize, Deserialize, PartialEq)]
        struct Input {
            #[serde(with = "scalar")]
            scalar: Rat,
            #[serde(with = "vector")]
            vector: Vec<Rat>,
            #[serde(with = "matrix")]
            matrix: Mat,
        }
        let input: Input = serde_json::from_str(r#"{"scalar":["-6","-9"],"vector":["1.25","9007199254740993"],"matrix":[["0",["1","3"]]]}"#).unwrap();
        assert_eq!(input.scalar, r("2/3"));
        assert_eq!(
            input,
            serde_json::from_str(&serde_json::to_string(&input).unwrap()).unwrap()
        );
        for invalid in [
            r#"{"scalar":0.1,"vector":[],"matrix":[]}"#,
            r#"{"scalar":[1,2],"vector":[],"matrix":[]}"#,
            r#"{"scalar":["1","0"],"vector":[],"matrix":[]}"#,
        ] {
            assert!(serde_json::from_str::<Input>(invalid).is_err());
        }
    }

    #[test]
    fn row_and_column_pivots_exact_residuals() {
        let a = m(&[&["0", "0", "2"], &["0", "3", "1"], &["5", "2", "0"]]);
        let b = vec![r("6"), r("9"), r("9")];
        let x = unique(&a, &b, 3).unwrap();
        assert_eq!(x, vec![r("1"), r("2"), r("3")]);
        for (row, rhs) in a.iter().zip(&b) {
            assert_eq!(row.iter().zip(&x).map(|(a, x)| a * x).sum::<Rat>(), *rhs);
        }
        let free = reduce(
            &m(&[&["0", "1", "2"], &["0", "2", "4"]]),
            &[r("3"), r("6")],
            3,
        );
        assert_eq!(free.pivots, vec![1]);
        assert!(!free.inconsistent);
        let kernel = [r("0"), r("-2"), r("1")];
        assert_eq!(
            free.rows[0]
                .iter()
                .zip(kernel)
                .map(|(a, x)| a * x)
                .sum::<Rat>(),
            Rat::zero()
        );
    }

    #[test]
    fn inconsistency_singularity_and_tiny_exact_pivot() {
        let a = m(&[&["1", "1"], &["2", "2"]]);
        assert!(reduce(&a, &[r("1"), r("3")], 2).inconsistent);
        assert!(unique(&a, &[r("1"), r("2")], 2).is_none());
        let tiny = r("1/1267650600228229401496703205376");
        let a = vec![vec![r("1"), r("1")], vec![r("1"), r("1") + &tiny]];
        let b = vec![r("3"), r("3") + r("2") * tiny];
        let reduced = reduce(&a, &b, 2);
        assert_eq!(reduced.rhs, vec![r("1"), r("2")]);
        assert_eq!(reduced.pivots, vec![0, 1]);
        assert!(reduced.peak_bits > 100);
        assert_eq!(rank(&a, 2), 2);
    }
}
