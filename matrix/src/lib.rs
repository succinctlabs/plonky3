//! Matrix library.

#![no_std]

extern crate alloc;

use alloc::boxed::Box;

pub mod dense;
pub mod mul;
pub mod sparse;
pub mod stack;

pub trait Matrix<T> {
    fn width(&self) -> usize;
    fn height(&self) -> usize;
}

/// A `Matrix` that supports randomly accessing particular rows or entries. Generally dense matrices
/// should implement this, but sparse matrices may not.
pub trait DenseMatrix<T>: Matrix<T> {
    type Row<'a>: IntoIterator<Item = &'a T>
    where
        Self: 'a,
        T: 'a;

    fn row<'a>(&'a self, r: usize) -> Self::Row<'a>;

    fn get(&self, r: usize, c: usize) -> T
    where
        T: Clone,
    {
        self.row(r)
            .into_iter()
            .nth(c)
            .expect("c out of range")
            .clone()
    }
}

impl<T> Matrix<T> for Box<dyn Matrix<T>> {
    fn width(&self) -> usize {
        self.as_ref().width()
    }

    fn height(&self) -> usize {
        self.as_ref().height()
    }

    // fn row(&self, r: usize) -> &[T] {
    //     self.as_ref().row(r)
    // }
}
