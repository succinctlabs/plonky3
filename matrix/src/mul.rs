use crate::dense::RowMajorMatrix;
use crate::sparse::CsrMatrix;
use crate::{DenseMatrix, Matrix};
use alloc::vec;
use p3_field::Field;

/// Compute `C = A * B`, where `A` in a CSR matrix and `B` is a dense matrix.
///
/// # Panics
/// Panics if dimensions of input matrices don't match.
pub fn mul_csr_dense_v2<F: Field, B: DenseMatrix<F>>(a: &CsrMatrix<F>, b: &B) -> RowMajorMatrix<F> {
    assert_eq!(a.width(), b.height(), "A, B dimensions don't match");
    let c_width = b.width();
    let c_height = a.height();
    let mut c = RowMajorMatrix::new(vec![F::ZERO; c_width * c_height], c_width);

    for a_row_idx in 0..a.height() {
        let c_row = c.row_mut(a_row_idx);
        for &(a_col_idx, a_val) in a.row(a_row_idx) {
            add_scaled_slice_in_place(c_row, b.row(a_col_idx).into_iter(), a_val);
        }
    }

    c
}

/// `x += y * s`, where `s` is a scalar.
fn add_scaled_slice_in_place<'a, F, Y>(x: &'a mut [F], y: Y, s: F)
where
    F: Field,
    Y: Iterator<Item = &'a F>,
{
    // TODO: Use PackedField
    x.iter_mut().zip(y).for_each(|(x_i, y_i)| *x_i += *y_i * s);
}
