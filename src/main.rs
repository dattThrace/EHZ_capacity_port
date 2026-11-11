use ehz_capacity::{ehz_from_facets, exact, geometry, Vertices};
use serde::Deserialize;
use serde_json::json;
use std::io::Read;
#[derive(Deserialize)]
#[serde(untagged)]
enum Input {
    Vertices {
        vertices: Vertices,
    },
    Facets {
        normals: Vertices,
        heights: Vec<f64>,
    },
}

// Inspect only the schema discriminator. In particular, exact scalar tokens are
// never first decoded into serde_json::Value's machine-number representation.
#[derive(Deserialize)]
struct InputSchema {
    schema: Option<String>,
}

fn emit(output: &serde_json::Value) -> Result<(), String> {
    println!(
        "{}",
        serde_json::to_string_pretty(output).map_err(|e| e.to_string())?
    );
    Ok(())
}

fn exact_error(error: String) -> Result<(), String> {
    emit(&json!({"result": {
        "schema": "ehz-exact-result-v1", "status": "invalid-input", "error": error
    }}))?;
    Err(error)
}

fn run() -> Result<(), String> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.iter().any(|a| a == "--help" || a == "-h") {
        println!(
            "ehz-capacity [--mode numerical|exact] [--geometry] [--max-steps N] [--timeout-ms N] [--certify-interval] [INPUT.json|-]\n\
             Numerical input: {{\"vertices\":[...]}} or {{\"normals\":[...],\"heights\":[...]}}.\n\
             Exact input: {{\"schema\":\"ehz-exact-h-v1\",\"normals\":[[\"1\",\"0\"],...],\"heights\":[\"1\",...]}}.\n\
             Exact schema selects exact mode automatically. Limits and --certify-interval apply to exact mode only.\n\
             --certify-interval adds a conservative certified interval to a partial positive witness.\n\
             Output: JSON. Indices: zero based. See docs/EXACT.md."
        );
        return Ok(());
    }
    let mut path = None;
    let mut geometry_only = false;
    let mut certify_interval = false;
    let mut mode = None;
    let mut options = exact::SearchOptions::default();
    let mut args = args.into_iter();
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--geometry" => geometry_only = true,
            "--certify-interval" => certify_interval = true,
            "--mode" => {
                let value = args.next().ok_or("--mode requires numerical or exact")?;
                if value != "numerical" && value != "exact" {
                    return Err(format!("unknown mode: {value}"));
                }
                if mode.replace(value).is_some() {
                    return Err("--mode can be supplied only once".into());
                }
            }
            "--max-steps" | "--timeout-ms" => {
                let value = args
                    .next()
                    .ok_or_else(|| format!("{arg} requires an integer"))?;
                let value = value
                    .parse::<u64>()
                    .map_err(|_| format!("{arg} requires a nonnegative integer"))?;
                let target = if arg == "--max-steps" {
                    &mut options.max_steps
                } else {
                    &mut options.timeout_ms
                };
                if target.replace(value).is_some() {
                    return Err(format!("{arg} can be supplied only once"));
                }
            }
            _ if arg.starts_with('-') && arg != "-" => {
                return Err(format!("unknown option: {arg}"));
            }
            _ => {
                if path.replace(arg).is_some() {
                    return Err("only one input file is supported".into());
                }
            }
        }
    }
    let mut input = String::new();
    if let Some(p) = path.filter(|p| p != "-") {
        input = std::fs::read_to_string(p).map_err(|e| e.to_string())?;
    } else {
        std::io::stdin()
            .read_to_string(&mut input)
            .map_err(|e| e.to_string())?;
    }
    let schema = match serde_json::from_str::<InputSchema>(&input) {
        Ok(header) => header.schema,
        Err(error) if mode.as_deref() == Some("exact") => {
            return exact_error(error.to_string());
        }
        Err(error) => return Err(error.to_string()),
    };
    if let Some(value) = schema.as_deref() {
        if value != "ehz-exact-h-v1" {
            return exact_error(format!("unsupported input schema: {value}"));
        }
        if mode.as_deref() == Some("numerical") {
            return exact_error("exact input schema conflicts with --mode numerical".into());
        }
    }
    if mode.as_deref() == Some("exact") || schema.is_some() {
        if geometry_only {
            return exact_error("--geometry is unavailable for exact H-input".into());
        }
        let exact_input = match serde_json::from_str::<exact::ExactInput>(&input) {
            Ok(value) => value,
            Err(error) => return exact_error(error.to_string()),
        };
        let mut result = exact::solve(&exact_input, &options);
        if certify_interval {
            exact::certify_interval(&mut result)?;
        }
        let result = json!(result);
        let invalid = result["status"] == "invalid-input";
        emit(&json!({"result": result}))?;
        return if invalid {
            Err("invalid exact input".into())
        } else {
            Ok(())
        };
    }
    if options.max_steps.is_some() || options.timeout_ms.is_some() || certify_interval {
        return Err("--max-steps, --timeout-ms, and --certify-interval require exact mode".into());
    }
    let output = match serde_json::from_str::<Input>(&input).map_err(|e| e.to_string())? {
        Input::Vertices { vertices } => {
            let g = geometry(&vertices)?;
            if geometry_only {
                json!(g)
            } else {
                let r = ehz_from_facets(&g.normals, &g.heights)?;
                json!({"geometry":g,"result":r})
            }
        }
        Input::Facets { normals, heights } => {
            if geometry_only {
                return Err("--geometry requires vertices".into());
            }
            json!({"result":ehz_from_facets(&normals,&heights)?})
        }
    };
    emit(&output)
}
fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error}");
        std::process::exit(1);
    }
}
