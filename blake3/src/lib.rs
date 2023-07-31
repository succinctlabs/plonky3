//! The BLAKE3 permutation, and hash functions built from it.

#![no_std]

extern crate alloc;

use alloc::vec::Vec;

use blake3::{hash_many_wrapper, MAX_SIMD_DEGREE_OR_2};
use itertools::Itertools;

/// The BLAKE3 permutation.
pub struct Blake3Hash;

const VSIZE: usize = MAX_SIMD_DEGREE_OR_2;

// impl CryptographicHasher<u8, [[u8; 32]; VSIZE]> for Blake3Hash {
impl Blake3Hash {
    pub fn hash_iter<I>(&self, input: I, input_len: usize) -> [[u8; 32]; VSIZE]
    where
        I: IntoIterator<Item = u8>,
    {
        let chunk_size = input_len / VSIZE;
        let input_vecs_vec: Vec<Vec<u8>> = input
            .into_iter()
            .chunks(chunk_size)
            .into_iter()
            .map(|chunk| chunk.collect_vec())
            .collect_vec();
        let input_slices_array = input_vecs_vec
            .iter()
            .map(|vec| vec.as_slice())
            .collect_vec()
            .try_into()
            .unwrap();

        let output = hash_many_wrapper(&input_slices_array);
        output
            .intern
            .chunks(chunk_size)
            .map(|chunk| chunk.try_into().unwrap())
            .collect_vec()
            .try_into()
            .unwrap()
    }

}
