//! The BLAKE3 permutation, and hash functions built from it.

#![no_std]

extern crate alloc;

use alloc::vec::Vec;
use blake3::{MAX_SIMD_DEGREE_OR_2, hash_many_wrapper};
use p3_symmetric::compression::PseudoCompressionFunction;

/// The BLAKE3 permutation.
pub struct Blake3Compression;

const COMPRESS_OUT_LEN: usize = 32 * MAX_SIMD_DEGREE_OR_2;

impl PseudoCompressionFunction<[u8; COMPRESS_OUT_LEN], 32> for Blake3Compression {
    fn compress(&self, input: [[u8; COMPRESS_OUT_LEN]; 32]) -> [u8; COMPRESS_OUT_LEN] {
        let input_flattened = input.into_iter().flatten().collect::<Vec<_>>();
        let mut input_arrays_vec: Vec<[u8; 1024]> = Vec::new();
        for idx in 0..MAX_SIMD_DEGREE_OR_2 {
            input_arrays_vec.push(input_flattened[idx * 1024..(idx + 1) * 1024].try_into().unwrap());
        }
        let input_arrays_vec_cloned = input_arrays_vec.clone();
        let input_slices_vec: Vec<&[u8]> = input_arrays_vec_cloned.iter().map(|a| a.as_ref()).collect();
        let input_resized: [&[u8]; MAX_SIMD_DEGREE_OR_2] = input_slices_vec.try_into().unwrap();

        let output = hash_many_wrapper(&input_resized);
        output.intern
    }
}
