use p3_field::AbstractField;
use p3_mds::MdsPermutation;
use p3_symmetric::Permutation;

/// A generic trait for the external linear layer of a Poseidon2 permutation.
pub trait ExternalLinearLayer<T: Clone, const WIDTH: usize>: Permutation<[T; WIDTH]> {}

#[derive(Clone, Default)]
pub struct CirculantMD4<Mds> {
    m4: Mds,
}

impl<AF: AbstractField, Mds, const WIDTH: usize> Permutation<[AF; WIDTH]> for CirculantMD4<Mds>
where
    Mds: MdsPermutation<AF, 4>,
{
    fn permute_mut(&self, input: &mut [AF; WIDTH]) {
        // First, apply Diag(M4, ... , M4)
        self.matmul_m4(input);

        // Apply the circulant matrix (
        let t4 = WIDTH / 4;
        let mut stored = [AF::zero(), AF::zero(), AF::zero(), AF::zero()];
        for l in 0..4 {
            stored[l] = input[l].clone();
            for j in 1..t4 {
                stored[l] += input[4 * j + l].clone();
            }
        }
        for i in 0..input.len() {
            input[i] += stored[i % 4].clone();
        }
    }
}

impl<Mds> CirculantMD4<Mds> {
    /// Applies the block diagonal matrix Diag(M4, ... , M4) to the input.
    ///
    /// This is equivalent to applying the matrix M4 to each 4 elements of the input.
    fn matmul_m4<T: Clone, const WIDTH: usize>(&self, input: &mut [T; WIDTH])
    where
        Mds: MdsPermutation<T, 4>,
    {
        input
            .chunks_exact_mut(4)
            .for_each(|x| self.m4.permute_mut(x.try_into().unwrap()));
    }
}
