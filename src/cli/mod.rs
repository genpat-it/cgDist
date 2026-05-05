// mod.rs - CLI module

pub mod args;
pub mod config;
pub mod merge;
pub mod validation;

// Re-export main types for convenience
pub use args::Args;
pub use config::Config;
pub use validation::{validate_args, ValidationResult};
