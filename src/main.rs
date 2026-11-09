use ehz_capacity::{ehz_from_facets, geometry, Vertices};
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
fn run() -> Result<(), String> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.iter().any(|a| a == "--help" || a == "-h") {
        println!("ehz-capacity [--geometry] [INPUT.json|-]\nInput: {{\"vertices\":[...]}} or {{\"normals\":[...],\"heights\":[...]}}. Output: JSON. Indices: zero based.");
        return Ok(());
    }
    let mut path = None;
    let mut geometry_only = false;
    for arg in args {
        if arg == "--geometry" {
            geometry_only = true;
        } else if arg.starts_with('-') && arg != "-" {
            return Err(format!("unknown option: {arg}"));
        } else if path.replace(arg).is_some() {
            return Err("only one input file is supported".into());
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
    println!(
        "{}",
        serde_json::to_string_pretty(&output).map_err(|e| e.to_string())?
    );
    Ok(())
}
fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error}");
        std::process::exit(1);
    }
}
