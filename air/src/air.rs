use core::ops::{Add, Mul, Sub};
use p3_field::{AbstractField, AbstractionOf, Field};
use p3_matrix::dense::RowMajorMatrix;
use p3_matrix::DenseMatrix;

pub trait Air<AB: AirBuilder>: Sync {
    fn eval(&self, builder: &mut AB);

    fn preprocessed_trace(&self) -> Option<RowMajorMatrix<AB::F>> {
        None
    }
}

pub trait AirBuilder: Sized {
    type F: Field;

    type Expr: AbstractionOf<Self::F>
        + Add<Self::Var, Output = Self::Expr>
        + Sub<Self::Var, Output = Self::Expr>
        + Mul<Self::Var, Output = Self::Expr>;

    type Var: 'static
        + Into<Self::Expr>
        + Copy
        + Add<Self::F, Output = Self::Expr>
        + Add<Self::Var, Output = Self::Expr>
        + Add<Self::Expr, Output = Self::Expr>
        + Sub<Self::F, Output = Self::Expr>
        + Sub<Self::Var, Output = Self::Expr>
        + Sub<Self::Expr, Output = Self::Expr>
        + Mul<Self::F, Output = Self::Expr>
        + Mul<Self::Var, Output = Self::Expr>
        + Mul<Self::Expr, Output = Self::Expr>;

    type M<'a>: DenseMatrix<Self::Var, Row<'a> = &'a [Self::Var]>
    where
        Self: 'a;

    fn main<'a>(&'a self) -> Self::M<'a>;

    fn is_first_row(&self) -> Self::Expr;
    fn is_last_row(&self) -> Self::Expr;
    fn is_transition(&self) -> Self::Expr {
        self.is_transition_window(2)
    }
    fn is_transition_window(&self, size: usize) -> Self::Expr;

    /// Returns a sub-builder whose constraints are enforced only when `condition` is nonzero.
    fn when<I: Into<Self::Expr>>(&mut self, condition: I) -> FilteredAirBuilder<Self> {
        FilteredAirBuilder {
            inner: self,
            condition: condition.into(),
        }
    }

    /// Returns a sub-builder whose constraints are enforced only when `x != y`.
    fn when_ne<I1: Into<Self::Expr>, I2: Into<Self::Expr>>(
        &mut self,
        x: I1,
        y: I2,
    ) -> FilteredAirBuilder<Self> {
        self.when(x.into() - y.into())
    }

    /// Returns a sub-builder whose constraints are enforced only on the first row.
    fn when_first_row(&mut self) -> FilteredAirBuilder<Self> {
        self.when(self.is_first_row())
    }

    /// Returns a sub-builder whose constraints are enforced only on the last row.
    fn when_last_row(&mut self) -> FilteredAirBuilder<Self> {
        self.when(self.is_last_row())
    }

    /// Returns a sub-builder whose constraints are enforced on all rows except the last.
    fn when_transition(&mut self) -> FilteredAirBuilder<Self> {
        self.when(self.is_transition())
    }

    /// Returns a sub-builder whose constraints are enforced on all rows except the last `size - 1`.
    fn when_transition_window(&mut self, size: usize) -> FilteredAirBuilder<Self> {
        self.when(self.is_transition_window(size))
    }

    fn assert_zero<I: Into<Self::Expr>>(&mut self, x: I);

    fn assert_one<I: Into<Self::Expr>>(&mut self, x: I) {
        self.assert_zero(x.into() - Self::Expr::ONE);
    }

    fn assert_eq<I1: Into<Self::Expr>, I2: Into<Self::Expr>>(&mut self, x: I1, y: I2) {
        self.assert_zero(x.into() - y.into());
    }

    /// Assert that `x` is a boolean, i.e. either 0 or 1.
    fn assert_bool<I: Into<Self::Expr>>(&mut self, x: I) {
        let x = x.into();
        self.assert_zero(x.clone() * (x - Self::Expr::ONE));
    }
}

pub trait PairBuilder: AirBuilder {
    fn preprocessed<'a>(&'a self) -> Self::M<'a>;
}

pub trait PermutationAirBuilder: AirBuilder {
    fn permutation<'a>(&'a self) -> Self::M<'a>;

    fn permutation_randomness(&self) -> &[Self::Expr];
}

pub struct FilteredAirBuilder<'a, AB: AirBuilder> {
    inner: &'a mut AB,
    condition: AB::Expr,
}

impl<'a, AB: AirBuilder> AirBuilder for FilteredAirBuilder<'a, AB> {
    type F = AB::F;
    type Expr = AB::Expr;
    type Var = AB::Var;
    type M<'b> = AB::M<'b> where Self: 'b;

    fn main<'b>(&'b self) -> Self::M<'b> {
        self.inner.main()
    }

    fn is_first_row(&self) -> Self::Expr {
        self.inner.is_first_row()
    }

    fn is_last_row(&self) -> Self::Expr {
        self.inner.is_last_row()
    }

    fn is_transition_window(&self, size: usize) -> Self::Expr {
        self.inner.is_transition_window(size)
    }

    fn assert_zero<I: Into<Self::Expr>>(&mut self, x: I) {
        self.inner.assert_zero(self.condition.clone() * x.into());
    }
}

#[cfg(test)]
mod tests {
    use crate::{Air, AirBuilder};
    use p3_matrix::DenseMatrix;

    struct FibonacciAir;

    impl<AB: AirBuilder> Air<AB> for FibonacciAir {
        fn eval<'a>(&'a self, builder: &'a mut AB) {
            let main: AB::M<'a> = builder.main();

            let x_0: AB::Var = main.row(0)[0];
            let x_1: AB::Var = main.row(1)[0];
            let x_2: AB::Var = main.row(2)[0];
            drop(main);

            builder.when_first_row().assert_zero(x_0);
            builder.when_first_row().assert_one(x_1);
            builder.when_transition().assert_eq(x_0 + x_1, x_2);
        }
    }
}
