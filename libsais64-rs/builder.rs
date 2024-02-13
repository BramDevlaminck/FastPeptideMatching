use std::env;
use std::error::Error;
use std::path::PathBuf;
use std::process::Command;

fn main() -> Result<(), Box<dyn Error>> {
    // compile the c library
    Command::new("cmake").args(["-DCMAKE_BUILD_TYPE=\"Release\"", "libsais", "-Blibsais"]).status()?;
    Command::new("make").args(["-C", "libsais"]).status()?;
    Command::new("mv").args(["liblibsais.a", "libsais"]).status()?;

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