use std::env;
use std::error::Error;
use std::fmt::{Display, Formatter};
use std::path::PathBuf;
use std::process::{Command, ExitStatus};

/// Custom error for compilation of the C library
#[derive(Debug)]
struct CompileError<'a> {
    command: &'a str,
    exit_code: Option<i32>,
}

impl<'a> Display for CompileError<'a> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let end_text = if let Some(code) = self.exit_code {
            format!("with exit code {}", code)
        } else {
            "without exit code".to_string()
        };
        let text = format!("Command with name `{}` failed {}", self.command, end_text);
        write!(f, "{}", text)
    }
}

impl<'a> Error for CompileError<'a> {}

/// Handles the exit statuses of the executed bash commands
///
/// # Arguments
/// * `name` - Name of the executed bash command
/// * `exit_states` - The exit status of the executed bash command
///
/// # Returns
///
/// Returns () if the exit status was success
/// 
/// # Errors
/// 
/// Returns a CompilationError if the command failed
fn exit_status_to_result(name: &str, exit_status: ExitStatus) -> Result<(), CompileError> {
    match exit_status.success() {
        true => Ok(()),
        false => Err(CompileError {
            command: name,
            exit_code: exit_status.code(),
        }),
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    // compile the c library
    Command::new("rm")
        .args(["libsais/CMakeCache.txt"])
        .status().unwrap_or_default(); // if removing fails, it is since the cmake cache did not exist, we just can ignore it
    exit_status_to_result(
        "cmake",
        Command::new("cmake")
            .args(["-DCMAKE_BUILD_TYPE=\"Release\"", "libsais", "-Blibsais"])
            .status()?,
    )?;
    exit_status_to_result(
        "make",
        Command::new("make").args(["-C", "libsais"]).status()?,
    )?;

    // link the c libsais library to rust
    println!("cargo:rustc-link-search=native=libsais64-rs/libsais");
    println!("cargo:rustc-link-lib=static=libsais");

    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header("libsais-wrapper.h")
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        // Finish the builder and generate the bindings.
        .generate()?;

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from(env::var("OUT_DIR")?);
    bindings.write_to_file(out_path.join("bindings.rs"))?;

    Ok(())
}
