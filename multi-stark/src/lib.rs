//! A minimal multivariate STARK framework.

#![no_std]

extern crate alloc;

mod config;
mod folder;
mod prover;

pub use config::*;
pub use folder::*;
pub use prover::*;
