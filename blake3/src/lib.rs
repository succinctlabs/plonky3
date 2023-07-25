//! The BLAKE3 permutation, and hash functions built from it.

#![no_std]

extern crate alloc;

use alloc::vec::Vec;
use blake3::{MAX_SIMD_DEGREE_OR_2, hash_many_wrapper};
use itertools::Itertools;
use p3_symmetric::hasher::CryptographicHasher;

/// The BLAKE3 permutation.
pub struct Blake3Hash;

const VSIZE: usize = MAX_SIMD_DEGREE_OR_2;

impl CryptographicHasher<u8, [[u8; 32]; VSIZE]> for Blake3Hash {
    fn hash_iter<I>(&self, input: I) -> [[u8; 32]; VSIZE]
    where
        I: IntoIterator<Item = u8>,
    {
        let input_vec = input.into_iter().collect::<Vec<_>>();
        let chunk_size = input_vec.len() / VSIZE;
        let input_chunks = input_vec.chunks(chunk_size);
        let input_chunks_array: [&[u8]; VSIZE] = input_chunks.map(|chunk| chunk.try_into().unwrap()).collect_vec().try_into().unwrap();

        let output = hash_many_wrapper(&input_chunks_array);
        output.intern.chunks(chunk_size).map(|chunk| chunk.try_into().unwrap()).collect_vec().try_into().unwrap()
    }

    fn hash_iter_slices<'a, I>(&self, _input: I) -> [[u8; 32]; VSIZE]
    where
        I: IntoIterator<Item = &'a [u8]>,
    {
        todo!()
    }
}
