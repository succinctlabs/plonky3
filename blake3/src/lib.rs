//! The blake3 hash function.

#![no_std]

extern crate alloc;

use alloc::vec::Vec;

use p3_symmetric::CryptographicHasher;

/// The blake3 hash function.
#[derive(Copy, Clone)]
pub struct Blake3;

impl CryptographicHasher<u8, [u8; 32]> for Blake3 {
    fn hash_iter<I>(&self, input: I) -> [u8; 32]
    where
        I: IntoIterator<Item = u8>,
    {
        let input = input.into_iter().collect::<Vec<_>>();
        self.hash_iter_slices([input.as_slice()])
    }

    fn hash_iter_slices<'a, I>(&self, input: I) -> [u8; 32]
    where
        I: IntoIterator<Item = &'a [u8]>,
    {
        let mut hasher = blake3::Hasher::new();
        for chunk in input.into_iter() {
            #[cfg(not(feature = "parallel"))]
            hasher.update(chunk);
            #[cfg(feature = "parallel")]
            hasher.update_rayon(chunk);
        }
        hasher.finalize().into()
    }
}

impl CryptographicHasher<u32, [u32; 8]> for Blake3 {
    fn hash_iter<I>(&self, input: I) -> [u32; 8]
    where
        I: IntoIterator<Item = u32>,
    {
        let mut input = input.into_iter().collect::<Vec<_>>();
        if input.len() <= blake3_zkvm::BLOCK_LEN {
            let size = input.len();
            input.resize(blake3_zkvm::BLOCK_LEN, 0u32);
            blake3_zkvm::hash_single_block(input.as_slice().try_into().unwrap(), size)
        } else {
            let ret = self.hash_iter_slices([input.as_slice()]);
            ret
        }
    }

    fn hash_iter_slices<'a, I>(&self, input: I) -> [u32; 8]
    where
        I: IntoIterator<Item = &'a [u32]>,
    {
        let mut zkvm_hasher = blake3_zkvm::Hasher::new();

        for chunk in input.into_iter() {
            zkvm_hasher.update(chunk);
        }
        let mut out: [u32; 8] = [0u32; 8];
        zkvm_hasher.finalize(&mut out);

        out
    }
}
