//! The Poseidon2 permutation.
//!
//! This implementation was based upon the following resources:
//! - https://github.com/HorizenLabs/poseidon2/blob/main/plain_implementations/src/poseidon2/poseidon2.rs
//! - https://eprint.iacr.org/2023/323.pdf

extern crate alloc;

mod babybear;
mod diffusion;
mod external;
mod goldilocks;
use alloc::vec::Vec;

pub use babybear::DiffusionMatrixBabybear;
pub use diffusion::DiffusionPermutation;
pub use external::*;
pub use goldilocks::DiffusionMatrixGoldilocks;
use p3_field::{AbstractField, PrimeField};
use p3_mds::m4::M4Mds;
use p3_symmetric::{CryptographicPermutation, Permutation};
use serde::{Deserialize, Serialize};
#[cfg(feature = "rand")]
use rand::distributions::Standard;
#[cfg(feature = "rand")]
use rand::prelude::Distribution;
#[cfg(feature = "rand")]
use rand::Rng;

const SUPPORTED_WIDTHS: [usize; 8] = [2, 3, 4, 8, 12, 16, 20, 24];

/// The Poseidon2 permutation.
#[derive(Clone)]
pub struct Poseidon2<F, Diffusion, const WIDTH: usize, const D: u64> {
    /// The number of external rounds.
    rounds_f: usize,

    /// The number of internal rounds.
    rounds_p: usize,

    /// The round constants.
    constants: Vec<[F; 16]>,

    /// The linear layer used in internal rounds (only needs diffusion property, not MDS).
    internal_linear_layer: Diffusion,

    /// The matrix `M4` used in the external linear layer.
    m4: M4Mds,
}

impl<F, Diffusion, const WIDTH: usize, const D: u64> Poseidon2<F, Diffusion, WIDTH, D>
where
    F: PrimeField,
{
    /// Create a new Poseidon2 configuration.
    pub fn new(
        rounds_f: usize,
        rounds_p: usize,
        constants: Vec<[F; WIDTH]>,
        internal_linear_layer: Diffusion,
    ) -> Self {
        assert!(SUPPORTED_WIDTHS.contains(&WIDTH));
        Self {
            rounds_f,
            rounds_p,
            constants,
            internal_linear_layer,
            m4: M4Mds,
        }
    }

    /// Create a new Poseidon2 configuration with random parameters.
    #[cfg(feature = "rand")]
    pub fn new_from_rng<R: Rng>(
        rounds_f: usize,
        rounds_p: usize,
        internal_mds: Diffusion,
        rng: &mut R,
    ) -> Self
    where
        Standard: Distribution<F>,
    {
        let mut constants = Vec::new();
        let rounds = rounds_f + rounds_p;
        for _ in 0..rounds {
            constants.push(rng.gen::<[F; WIDTH]>());
        }

        Self {
            rounds_f,
            rounds_p,
            constants,
            internal_linear_layer: internal_mds,
            m4: M4Mds,
        }
    }

    #[inline]
    fn add_rc<AF>(&self, state: &mut [AF; 16], rc: &[AF::F; 16])
    where
        AF: AbstractField<F = F>,
    {
        state
            .iter_mut()
            .zip(rc)
            .for_each(|(a, b)| *a += AF::from_f(*b));
    }

    #[inline]
    fn sbox_p<AF>(&self, input: &AF) -> AF
    where
        AF: AbstractField<F = F>,
    {
        input.exp_const_u64::<D>()
    }

    #[inline]
    fn sbox<AF>(&self, state: &mut [AF; 16])
    where
        AF: AbstractField<F = F>,
    {
        state.iter_mut().for_each(|a| *a = self.sbox_p(a));
    }

    #[inline]
    fn external_linear_permute_mut<AF>(&self, state: &mut [AF; WIDTH])
    where
        AF: AbstractField<F = F>,
    {
        matmul_external_mut(state, &self.m4);
    }
}

impl<AF, Diffusion, const WIDTH: usize, const D: u64> Permutation<[AF; WIDTH]>
    for Poseidon2<AF::F, Diffusion, WIDTH, D>
where
    AF: AbstractField,
    AF::F: PrimeField,
    Diffusion: DiffusionPermutation<AF, WIDTH>,
{
    fn permute_mut(&self, state: &mut [AF; 16]) {
        // The initial linear layer.
        self.external_linear_permute_mut(state);

        // The first half of the external rounds.
        let rounds = self.rounds_f + self.rounds_p;
        let rounds_f_beggining = self.rounds_f / 2;
        for r in 0..rounds_f_beggining {
            self.add_rc(state, &self.constants[r]);
            self.sbox(state);
            self.external_linear_permute_mut(state);
        }

        // The internal rounds.
        let p_end = rounds_f_beggining + self.rounds_p;
        for r in self.rounds_f..p_end {
            state[0] += AF::from_f(self.constants[r][0]);
            state[0] = self.sbox_p(&state[0]);
            self.internal_linear_layer.permute_mut(state);
        }

        // The second half of the external rounds.
        for r in p_end..rounds {
            self.add_rc(state, &self.constants[r]);
            self.sbox(state);
            self.external_linear_permute_mut(state);
        }
    }
}

impl<AF, Diffusion, const WIDTH: usize, const D: u64> CryptographicPermutation<[AF; WIDTH]>
    for Poseidon2<AF::F, Diffusion, WIDTH, D>
where
    AF: AbstractField,
    AF::F: PrimeField,
    Diffusion: DiffusionPermutation<AF, WIDTH>,
{
}
