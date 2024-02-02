use p3_field::AbstractField;
use p3_mds::m4::M4Mds;
use p3_symmetric::Permutation;

pub fn matmul_external_mut<AF: AbstractField, const WIDTH: usize>(
    input: &mut [AF; WIDTH],
    m4: &M4Mds,
) {
    match WIDTH {
        2 => {
            // Matrix circ(2, 1)
            let mut sum = input[0].clone();
            sum += input[1].clone();
            input[0] += sum.clone();
            input[1] += sum.clone();
        }
        3 => {
            // Matrix circ(2, 1, 1)
            let mut sum = input[0].clone();
            sum += input[1].clone();
            sum += input[2].clone();
            input[0] += sum.clone();
            input[1] += sum.clone();
            input[2] += sum.clone();
        }
        4 => {
            // Applying cheap 4x4 MDS matrix to each 4-element part of the state
            matmul_m4(input, m4)
        }
        8 | 12 | 16 | 20 | 24 => {
            // First, apply Diag(M4, ... , M4)
            matmul_m4(input, m4);

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
        _ => panic!("Unsupported width"),
    }
}

#[inline]
fn matmul_m4<AF: AbstractField, const WIDTH: usize>(input: &mut [AF; WIDTH], m4: &M4Mds) {
    input
        .chunks_exact_mut(4)
        .for_each(|x| m4.permute_mut(x.try_into().unwrap()));
}
