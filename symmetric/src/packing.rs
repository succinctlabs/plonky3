use core::{mem, slice};

use p3_field::PackedField;

/// # Safety
/// - `WIDTH` is assumed to be a power of 2.
/// - If `P` implements `PackedField` then `P` must be castable to/from `[P::Value; P::WIDTH]`
///   without UB.
pub unsafe trait PackedValue: 'static + Copy + From<Self::Value> + Send + Sync {
    type Value: 'static + Copy + Send + Sync;

    const WIDTH: usize;

    fn from_slice(slice: &[Self::Value]) -> &Self;
    fn from_slice_mut(slice: &mut [Self::Value]) -> &mut Self;

    /// Similar to `core:array::from_fn`.
    fn from_fn<F>(f: F) -> Self
    where
        F: FnMut(usize) -> Self::Value;

    fn as_slice(&self) -> &[Self::Value];
    fn as_slice_mut(&mut self) -> &mut [Self::Value];

    fn pack_slice(buf: &[Self::Value]) -> &[Self] {
        // Sources vary, but this should be true on all platforms we care about.
        // This should be a const assert, but trait methods can't access `Self` in a const context,
        // even with inner struct instantiation. So we will trust LLVM to optimize this out.
        assert!(mem::align_of::<Self>() <= mem::align_of::<Self::Value>());
        assert!(
            buf.len() % Self::WIDTH == 0,
            "Slice length (got {}) must be a multiple of packed field width ({}).",
            buf.len(),
            Self::WIDTH
        );
        let buf_ptr = buf.as_ptr().cast::<Self>();
        let n = buf.len() / Self::WIDTH;
        unsafe { slice::from_raw_parts(buf_ptr, n) }
    }

    fn pack_slice_with_suffix(buf: &[Self::Value]) -> (&[Self], &[Self::Value]) {
        let (packed, suffix) = buf.split_at(buf.len() - buf.len() % Self::WIDTH);
        (Self::pack_slice(packed), suffix)
    }

    fn pack_slice_mut(buf: &mut [Self::Value]) -> &mut [Self] {
        assert!(mem::align_of::<Self>() <= mem::align_of::<Self::Value>());
        assert!(
            buf.len() % Self::WIDTH == 0,
            "Slice length (got {}) must be a multiple of packed field width ({}).",
            buf.len(),
            Self::WIDTH
        );
        let buf_ptr = buf.as_mut_ptr().cast::<Self>();
        let n = buf.len() / Self::WIDTH;
        unsafe { slice::from_raw_parts_mut(buf_ptr, n) }
    }

    fn unpack_slice(buf: &[Self]) -> &[Self::Value] {
        assert!(mem::align_of::<Self>() >= mem::align_of::<Self::Value>());
        let buf_ptr = buf.as_ptr().cast::<Self::Value>();
        let n = buf.len() * Self::WIDTH;
        unsafe { slice::from_raw_parts(buf_ptr, n) }
    }
}

unsafe impl<P: PackedField> PackedValue for P {
    type Value = P::Scalar;

    const WIDTH: usize = P::WIDTH;

    fn from_slice(slice: &[Self::Value]) -> &Self {
        P::from_slice(slice)
    }

    fn from_slice_mut(slice: &mut [Self::Value]) -> &mut Self {
        P::from_slice_mut(slice)
    }

    fn from_fn<F>(f: F) -> Self
    where
        F: FnMut(usize) -> Self::Value,
    {
        P::from_fn(f)
    }

    fn as_slice(&self) -> &[Self::Value] {
        P::as_slice(self)
    }

    fn as_slice_mut(&mut self) -> &mut [Self::Value] {
        P::as_slice_mut(self)
    }
}
