use p3_matrix::{DenseMatrix, Matrix};

#[derive(Copy, Clone)]
pub struct TwoRowMatrixView<'a, T> {
    pub local: &'a [T],
    pub next: &'a [T],
}

impl<'a, T> TwoRowMatrixView<'a, T> {
    pub fn new(local: &'a [T], next: &'a [T]) -> Self {
        Self { local, next }
    }
}

impl<'a, T> Matrix<T> for TwoRowMatrixView<'a, T> {
    fn width(&self) -> usize {
        self.local.len()
    }

    fn height(&self) -> usize {
        2
    }
}

impl<T> DenseMatrix<T> for TwoRowMatrixView<'_, T> {
    type Row<'a> = &'a [T] where Self: 'a, T: 'a;

    fn row<'a>(&'a self, r: usize) -> &'a [T] {
        match r {
            0 => self.local,
            1 => self.next,
            _ => panic!("Only two rows available"),
        }
    }
}
