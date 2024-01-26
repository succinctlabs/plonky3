use core::{
    marker::PhantomData,
    ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::{AbstractExtensionField, AbstractField, ExtensionField, Field};

/// The restriction of scalars from a field `EF` to a subfield `F`.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Default, Hash)]
pub struct Res<F: Field, EF: AbstractExtensionField<F>>(EF, PhantomData<F>);

impl<F: Field, EF: AbstractExtensionField<F>> Res<F, EF> {}

impl<F: Field, EF: AbstractExtensionField<F>> AbstractField for Res<F, EF> {
    type F = F;

    fn from_f(f: F) -> Self {
        EF::from_base(f).into()
    }

    fn zero() -> Self {
        EF::zero().into()
    }

    fn one() -> Self {
        EF::one().into()
    }

    fn two() -> Self {
        EF::two().into()
    }

    fn from_bool(b: bool) -> Self {
        EF::from_bool(b).into()
    }

    fn from_canonical_u8(n: u8) -> Self {
        EF::from_canonical_u8(n).into()
    }

    fn from_canonical_u16(n: u16) -> Self {
        EF::from_canonical_u16(n).into()
    }

    fn from_canonical_u32(n: u32) -> Self {
        EF::from_canonical_u32(n).into()
    }

    fn from_canonical_u64(n: u64) -> Self {
        EF::from_canonical_u64(n).into()
    }

    fn from_canonical_usize(n: usize) -> Self {
        EF::from_canonical_usize(n).into()
    }

    fn from_wrapped_u32(n: u32) -> Self {
        EF::from_wrapped_u32(n).into()
    }

    fn from_wrapped_u64(n: u64) -> Self {
        EF::from_wrapped_u64(n).into()
    }

    fn neg_one() -> Self {
        EF::neg_one().into()
    }

    fn generator() -> Self {
        EF::generator().into()
    }
}

impl<F: Field, EF: AbstractExtensionField<F>> From<EF> for Res<F, EF> {
    fn from(ef: EF) -> Self {
        Res(ef, PhantomData)
    }
}

impl<F: Field, EF: AbstractExtensionField<F>> Add for Res<F, EF> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        (self.0 + other.0).into()
    }
}

impl<F: Field, EF: AbstractExtensionField<F>> AddAssign for Res<F, EF> {
    fn add_assign(&mut self, other: Self) {
        self.0.add_assign(other.0);
    }
}

impl<F: Field, EF: AbstractExtensionField<F>> Mul for Res<F, EF> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        (self.0 * other.0).into()
    }
}

impl<F: Field, EF: AbstractExtensionField<F>> MulAssign for Res<F, EF> {
    fn mul_assign(&mut self, other: Self) {
        self.0.mul_assign(other.0);
    }
}

impl<F: Field, EF: AbstractExtensionField<F>> Sub for Res<F, EF> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        (self.0 - other.0).into()
    }
}

impl<F: Field, EF: AbstractExtensionField<F>> Neg for Res<F, EF> {
    type Output = Self;

    fn neg(self) -> Self {
        (-self.0).into()
    }
}

impl<F: Field, EF: AbstractExtensionField<F>> SubAssign for Res<F, EF> {
    fn sub_assign(&mut self, other: Self) {
        self.0.sub_assign(other.0);
    }
}

impl<F: Field, EF: ExtensionField<F>> Div for Res<F, EF> {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        (self.0 / other.0).into()
    }
}

impl<F: Field, EF: AbstractExtensionField<F>> core::iter::Product for Res<F, EF> {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<F: Field, EF: AbstractExtensionField<F>> core::iter::Sum for Res<F, EF> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}
