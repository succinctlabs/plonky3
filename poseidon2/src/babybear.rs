//! Diffusion matrices for Babybear16 and Babybear24.
//!
//! Reference: https://github.com/HorizenLabs/poseidon2/blob/main/plain_implementations/src/poseidon2/poseidon2_instance_babybear.rs

use p3_baby_bear::BabyBear;
// use p3_baby_bear::IN_HASH;
use p3_field::{AbstractField, PrimeField32};
use p3_symmetric::Permutation;

use succinct_zkvm::{io, unconstrained};

use crate::diffusion::matmul_internal;
use crate::DiffusionPermutation;

pub const MATRIX_DIAG_16_BABYBEAR: [u64; 16] = [
    0x0a632d94, 0x6db657b7, 0x56fbdc9e, 0x052b3d8a, 0x33745201, 0x5c03108c, 0x0beba37b, 0x258c2e8b,
    0x12029f39, 0x694909ce, 0x6d231724, 0x21c3b222, 0x3c0904a5, 0x01d6acda, 0x27705c83, 0x5231c802,
];

pub const MATRIX_DIAG_24_BABYBEAR: [u64; 24] = [
    0x409133f0, 0x1667a8a1, 0x06a6c7b6, 0x6f53160e, 0x273b11d1, 0x03176c5d, 0x72f9bbf9, 0x73ceba91,
    0x5cdef81d, 0x01393285, 0x46daee06, 0x065d7ba6, 0x52d72d6f, 0x05dd05e0, 0x3bab4b63, 0x6ada3842,
    0x2fc5fbec, 0x770d61b0, 0x5715aae9, 0x03ef0e90, 0x75b6c770, 0x242adf5f, 0x00d0ca4c, 0x36c0e388,
];

#[derive(Debug, Clone, Default)]
pub struct DiffusionMatrixBabybear;

impl<AF: PrimeField32> Permutation<[AF; 16]> for DiffusionMatrixBabybear {
    fn permute_mut(&self, state: &mut [AF; 16]) {
        // let mut in_hash = IN_HASH.lock().unwrap();
        // *in_hash = true;
        // drop(in_hash);
        // println!("cycle-tracker-start: permute_mut matmul_internal");

        unconstrained! {
            let mut new_state: [AF; 16] = [AF::default(); 16];
            new_state.copy_from_slice(state);
            matmul_internal::<AF, 16>(&mut new_state, MATRIX_DIAG_16_BABYBEAR);
            let bytes = state.map(|x| x.as_canonical_u32().to_le_bytes());
            let mut flat_bytes = Vec::new();
            for i in 0..16 {
                flat_bytes.extend_from_slice(&bytes[i]);
            }
            io::hint_slice(&flat_bytes);
        }

        let mut bytes: [u8; 64] = [0; 64];
        io::read_hint_slice(&mut bytes);
        let ret = bytes.chunks(4).map(|chunk| AF::from_canonical_u32(u32::from_le_bytes(chunk.try_into().unwrap()))).collect::<Vec<AF>>();
        for i in 0..16 {
            state[i] = ret[i];
        }

        // println!("cycle-tracker-end: permute_mut matmul_internal");
        // let mut in_hash = IN_HASH.lock().unwrap();
        // *in_hash = false;
        // drop(in_hash);
    }
}

impl<AF: PrimeField32> DiffusionPermutation<AF, 16> for DiffusionMatrixBabybear {}

impl<AF: AbstractField<F = BabyBear>> Permutation<[AF; 24]> for DiffusionMatrixBabybear {
    fn permute_mut(&self, state: &mut [AF; 24]) {
        matmul_internal::<AF, 24>(state, MATRIX_DIAG_24_BABYBEAR);
    }
}

impl DiffusionPermutation<BabyBear, 24> for DiffusionMatrixBabybear {}
