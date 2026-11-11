use ehz_capacity::exact::arithmetic::{parse_rat, Rat};
use ehz_capacity::exact::{check_witness, omega, solve, ExactInput, ExactResult, SearchOptions};
use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::io::Write;
use std::process::{Command, Stdio};

fn square() -> ExactInput {
    serde_json::from_str(include_str!("../examples/exact/square.json")).unwrap()
}

fn fraction(text: &str) -> Rat {
    parse_rat(text).unwrap()
}

#[test]
fn square_requires_full_support_and_preserves_order_conventions() {
    let input = square();
    let result = solve(&input, &SearchOptions::default());
    assert_eq!(result.status, "exact-complete");
    assert!(result.completion.complete);
    assert_eq!(result.capacity, Some(fraction("4")));
    let witness = result.witness.unwrap();
    assert_eq!(witness.order.len(), 4);
    assert_eq!(witness.objective, fraction("1/8"));
    check_witness(&input, &witness).unwrap();

    let mut cyclic = witness.clone();
    cyclic.order.rotate_left(1);
    cyclic.weights.rotate_left(1);
    check_witness(&input, &cyclic).unwrap();
    let mut reverse_score = Rat::zero();
    cyclic.order.reverse();
    cyclic.weights.reverse();
    for i in 0..cyclic.order.len() {
        for j in 0..i {
            reverse_score += &cyclic.weights[i]
                * &cyclic.weights[j]
                * omega(
                    &input.normals[cyclic.order[i]],
                    &input.normals[cyclic.order[j]],
                );
        }
    }
    assert_eq!(reverse_score, -witness.objective);
    assert!(check_witness(&input, &cyclic).is_err());
}

#[test]
fn certificate_roundtrip_preserves_original_scaled_row_witness() {
    let input: ExactInput =
        serde_json::from_str(include_str!("../examples/exact/large_row_scales.json")).unwrap();
    let result = solve(&input, &SearchOptions::default());
    assert_eq!(result.status, "exact-complete");
    assert_eq!(result.capacity, Some(fraction("4")));
    let printed = serde_json::to_string(&result).unwrap();
    let recovered: ExactResult = serde_json::from_str(&printed).unwrap();
    assert_eq!(recovered.input, input);
    assert_eq!(recovered.witness, result.witness);
    let witness = recovered.witness.unwrap();
    check_witness(&input, &witness).unwrap();
    let mut mutated = witness.clone();
    mutated.weights[0] += fraction("1/100000000000000000000000000000000000000000");
    assert!(check_witness(&input, &mutated).is_err());
    mutated = witness.clone();
    mutated.objective += fraction("1/100000000000000000000000000000000000000000");
    assert!(check_witness(&input, &mutated).is_err());
}

#[test]
fn certified_redundancy_retains_the_same_square() {
    let input: ExactInput =
        serde_json::from_str(include_str!("../examples/exact/redundant_square.json")).unwrap();
    let result = solve(&input, &SearchOptions::default());
    assert_eq!(result.capacity, Some(fraction("4")));
    assert_eq!(result.status, "exact-complete");
    let domain = result.domain.unwrap();
    assert_eq!(domain.retained.len(), 4);
    assert_eq!(domain.rows.len(), input.normals.len());
    assert!(!domain.removals.is_empty());
    check_witness(&input, &result.witness.unwrap()).unwrap();
}

#[test]
fn rational_perturbations_survive_beyond_binary64() {
    let mut input = square();
    for exponent in [20, 60, 120] {
        let epsilon = Rat::new(BigInt::one(), BigInt::one() << exponent);
        for sign in [-1, 0, 1] {
            let delta = &epsilon * Rat::from_integer(sign.into());
            input.heights[0] = Rat::one() + &delta;
            let result = solve(&input, &SearchOptions::default());
            assert_eq!(result.status, "exact-complete");
            assert_eq!(result.capacity, Some(fraction("4") + fraction("2") * delta));
        }
    }
}

#[test]
fn interruption_never_promotes_a_partial_incumbent() {
    let input = square();
    for options in [
        SearchOptions {
            max_steps: Some(0),
            timeout_ms: None,
        },
        SearchOptions {
            max_steps: None,
            timeout_ms: Some(0),
        },
    ] {
        let result = solve(&input, &options);
        assert_eq!(result.status, "incomplete");
        assert!(!result.completion.complete);
        assert!(result.capacity.is_none());
        assert!(result.witness.is_none());
    }
    let complete = solve(&input, &SearchOptions::default());
    let result = solve(
        &input,
        &SearchOptions {
            max_steps: Some(complete.stats.steps - 1),
            timeout_ms: None,
        },
    );
    assert_eq!(result.status, "exact-witness/incomplete");
    assert!(!result.completion.complete);
    assert!(result.capacity.is_none());
    check_witness(&input, result.witness.as_ref().unwrap()).unwrap();
}

#[test]
fn exact_cli_rejects_machine_numbers_as_structured_invalid_input() {
    let mut input = serde_json::to_value(square()).unwrap();
    input["heights"][0] = serde_json::json!(1.0);
    let mut child = Command::new(env!("CARGO_BIN_EXE_ehz-capacity"))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(serde_json::to_string(&input).unwrap().as_bytes())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(!output.status.success());
    let result: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(result["result"]["status"], "invalid-input");
}
