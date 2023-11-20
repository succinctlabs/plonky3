//! An implementation of split-radix FFT for `Mersenne31` and
//! `Mersenne31Complex`.
//!
//! This code is essentially a port of the excellent `djbfft`
//! library[1] to Rust in for `Mersenne31Complex`. The main
//! "innovation" of the port is to introduce delayed modular reduction
//! wherever possible.
//!
//! NB: `djbfft` has not been updated in a long time, and is
//! specifically optimised for older processors. This port has not
//! (yet) made any special effort to adapt the instruction scheduling
//! to a modern processor, but we did of course verify that
//! performance is sufficient for our purposes. In some places where
//! choices had to be made, we have deliberately chosen clarity over
//! literal translation of the code, even when this might possibly
//! have compromised performance. Future performance-oriented
//! travellers to these lands should therefore not take the
//! instruction layout as gospel.
//!
//! [1] https://cr.yp.to/djbfft.html

// Implements forward FFT using DIF and split-radix conjugate pair FFT.
//
// Almost everything in djbfft has a name built of single
// letters/digits that decode to describe what it does. This is
// essentially the result of manually specialising everything due to
// the lack of language support for "generic programming" in C.
// Anyway, the letter/digit decoding is:
//
// '4': 4-byte (single word) data
// '8': 8-byte (double word) data
// 'c': complex DIF
// 'r': real DIF
// 'd': roots of unity
// 'i': include file containing prototypes and macros
// 's': scaling
// 'u': complex DIT
// 'v': real DIT
//
// Sometimes a logical group of functionality is split between
// different files, with an incrementing number at the end used to
// distinguish the files. I'm just guessing, but I presume this was to
// allow independent compilation of object files, so changes in one
// wouldn't require the slow recompilation of everything.
//
// So, for example, the DIF split-radix complex FFT for 4-byte
// real/imaginary parts is implemented in the files `4c0.c`, `4c1.c`,
// `4c2.c`, `4c3.c`, `4c4.c`, and `4c5.c`, with all the small
// specialisations in `4c0.c`, new code for bigger specialisations in
// `4c1.c`, and then brief specialisations for sizes 1024, 2048, 4096
// and 8192 in the files `4c2.c`...`4c5.c` respectively.

use core::fmt;

use crate::split_radix::complex_backward::*;
use crate::split_radix::complex_forward::*;

// TODO: These are only pub for benches at the moment...
pub mod complex_backward;
pub mod complex_forward;
pub mod real;
pub(crate) mod roots;

// TODO: Not sure that forcing alignment is necessary...
//#[repr(align(8))]
//pub struct Real(i64);
type Real = i64;

#[derive(Copy, Clone, Default)]
#[repr(align(8))]
pub struct Complex {
    pub re: Real,
    pub im: Real,
}

pub fn forward_fft<const N: usize>(u: &mut [Complex]) {
    match N {
        8 => complex_forward_8(u),
        16 => complex_forward_16(u),
        32 => complex_forward_32(u),
        64 => complex_forward_64(u),
        128 => complex_forward_128(u),
        256 => complex_forward_256(u),
        512 => complex_forward_512(u),
        1024 => complex_forward_1024(u),
        2048 => complex_forward_2048(u),
        4096 => complex_forward_4096(u),
        8192 => complex_forward_8192(u),
        _ => panic!("unsupported size"),
    }
}

pub fn backward_fft<const N: usize>(v: &mut [Complex]) {
    match N {
        8 => complex_backward_8(v),
        16 => complex_backward_16(v),
        32 => complex_backward_32(v),
        64 => complex_backward_64(v),
        128 => complex_backward_128(v),
        256 => complex_backward_256(v),
        512 => complex_backward_512(v),
        1024 => complex_backward_1024(v),
        2048 => complex_backward_2048(v),
        4096 => complex_backward_4096(v),
        8192 => complex_backward_8192(v),
        _ => panic!("unsupported size"),
    }
}

const P: Real = (1 << 31) - 1;

impl Complex {
    pub const ZERO: Self = Self::new(0i64, 0i64);
    //const ONE: Self = Self::new(1i64, 0i64);
    //const J: Self = Self::new(0i64, 1i64);

    pub const fn new(re: Real, im: Real) -> Self {
        Self { re, im }
    }

    pub const fn from_pair((re, im): (Real, Real)) -> Self {
        Self::new(re, im)
    }
}

impl fmt::Debug for Complex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} + {}*J", self.re, self.im)
    }
}

pub(crate) const fn into<const N: usize>(u: [(Real, Real); N]) -> [Complex; N] {
    let mut v = [Complex::ZERO; N];
    let mut i = 0;
    loop {
        v[i] = Complex::from_pair(u[i]);
        i += 1;
        if i == N {
            break;
        }
    }
    v
}

#[inline(always)]
pub(crate) fn reduce_2p(x: Real) -> Real {
    debug_assert!(x >= 0 && x <= 2 * P, "x = {} out of range", x);

    // let msb = x & (1 << 31);
    // let msb_reduced = msb >> 31;
    // (x ^ msb) + msb_reduced

    // This method seems to make the FFT bench 5--9% faster than the
    // one above
    let x = x as u64;
    x.min(x.overflowing_sub(P as u64).0) as i64
}

/*
/// Given 0 <= x <= P, return x' such that
/// -(p - 1)/2 <= x' <= (p + 1)/2 and x = x' (mod p).
#[inline]
fn shift(x: i64) -> i64 {
    debug_assert!(x >= 0 && x <= P);

// As a u32, the top two bits of x are b = 00 or 01. In the first case,
// x < P/2 and so can be left alone. In the second case x > P/2 and we
// want to subtract P = 2^31 - 1 to get x - 2^31 + 1. As an i32, -2^31
// corresponds to the most significant bit, whence the "| (b<<31)", then
// we add 1.

    let x = x as u64;
    let b = (x >> 30) as i32;
    let y = (x as i32 | (b << 31)) + b;
    y as i64 // sign extend...
}
*/

/*
#[inline(always)]
fn reduce_4p_real(x: Real) -> Real {
    debug_assert!(x >= -3 * P && x <= 4 * P);
    const MASK: u64 = (1 << 31) - 1;
    let y = (x + 3 * P) as u64; // make unsigned so the right shift doesn't sign extend
    let (lo, hi) = ((y & MASK) as i64, ((y >> 31) & MASK) as i64);
    lo + hi
}

#[inline(always)]
pub(crate) fn reduce_4p(z: &mut Complex) {
    z.re = reduce_4p_real(z.re);
    z.im = reduce_4p_real(z.im);
}
*/

/// Given a Real x in the range |x| <= 2 * (P^2 - 1), return
/// a Real x' in the range 0 <= x' <= P such that x' = x (mod P).
#[inline(always)]
pub(crate) fn reduce_2p_sqr(x: Real) -> Real {
    const ABSMAX_ELT: i64 = 2 * (P * P - 1);
    debug_assert!(x <= ABSMAX_ELT && x >= -ABSMAX_ELT);

    const MASK: u64 = (1 << 31) - 1;
    let y = x as u64; // make unsigned so the right shift doesn't sign extend
    let (lo, hi) = ((y & MASK) as i64, ((y >> 31) & MASK) as i64);

    // Here we prove the correctness of the line "lo + hi + (x >> 62)".
    //
    // FIXME: Note that for one particular value (3 << 62) we "underflow"
    // and pass -1 to `normalise_u32()`. This is a bug.
    //
    // The input x is an i64. If we label its bits from b[0] to b[63]
    // and use the relation 2^31 = 1 (mod p), then
    //
    //   x = -b[63] * 2^63 + \sum_{i=0}^62 b[i] * 2^i  // 2's compl repr
    //     = -b[63] * 2^63 + b[62] * 2^62
    //          + \sum_{i=0}^30 b[i] * 2^i + \sum_{i=31}^61 b[i] * 2^{i-31}
    //     = -2 * b[63] + b[62] + b[0..=30] + b[31..=61]
    //     = -2 * b[63] + b[62] + lo + hi
    //
    // Note:
    //  1) 0 <= lo, hi <= 2^31 - 1 == P
    //  2) The possible values of -2 * b[63] + b[62] are:
    //        0 if b[63] = b[62] = 0
    //        1 if b[63] = 0 and b[62] = 1
    //       -2 if b[63] = 1 and b[62] = 0
    //       -1 if b[63] = b[62] = 1
    //     This is a signed 2-bit number in 2's complement, so can be
    //     obtained as the arithmetic right shift of the top two bits,
    //     whence (x >> 62) = -2 * b[63] + b[62].
    //  3) Regarding overflow and underflow, consider each possible value
    //     of -2 * b[63] + b[62]:
    //       a.  0 ==> Then the result is lo + hi which is fine.
    //       b.  1 ==> Would be > 2P if lo = hi = P, however that corresponds
    //              to the bit pattern b[63] = 0 and b[i] = 1 for 0 <= i < 63
    //              which is 2^63 - 1 > 2 * P * (P - 1), hence can't happen.
    //       c. -2 ==> Would be < 0 if lo + hi = 0 or 1. These possibilites
    //              correspond to the input x being:
    //                 0x8000_0000_0000_0000   lo = hi = 0
    //                 0x8000_0000_0000_0001   lo = 1, hi = 0
    //                 0x8000_0000_8000_0000   lo = 0, hi = 1
    //              all of which are < -2 * P * (P - 1), hence can't happen.
    //       d. -1 ==> Will be < 0 if lo = hi = 0. This corresponds to
    //              x = 0xC000_0000_0000_0000 > -2 * P * (P - 1), which
    //              is a valid input (specifically it is = -1 (mod p)).
    //              *** HENCE THIS IS A BUG!! ***
    //

    // 0 <= lo + hi <= 2P
    let fold = reduce_2p(lo + hi);
    // 0 <= fold < P (or == P if lo + hi = 2P)

    // x must be signed for the SAR here
    let x_mod_p = fold + (x >> 62);
    // -1 <= x_mod_p <= P

    // TODO: Work out whether it's okay to allow the bigger range and
    // avoid the second u32 reduction. One problem is that a0 and a1
    // in transform(...) could be -2 on output, so we'd need an
    // alternative to normalise_u32.

    //x_mod_p

    // P - 1 <= x_mod_p <= 2P
    reduce_2p(x_mod_p + P)
    // 0 <= x_mod_p <= P
}

#[inline(always)]
pub(crate) fn reduce_complex(z: &mut Complex) {
    z.re = reduce_2p_sqr(z.re);
    z.im = reduce_2p_sqr(z.im);
}

#[inline]
pub(crate) fn reduce_all(zs: &mut [Complex]) {
    zs.iter_mut().for_each(reduce_complex);
}

impl PartialEq for Complex {
    fn eq(&self, other: &Self) -> bool {
        reduce_2p_sqr(self.re) == reduce_2p_sqr(other.re)
            && reduce_2p_sqr(self.im) == reduce_2p_sqr(other.im)
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;
    use core::iter::repeat_with;
    use core::ops::{Add, Mul};

    use rand::{thread_rng, Rng};

    use crate::split_radix::complex_backward::{
        complex_backward_128, complex_backward_16, complex_backward_256, complex_backward_32,
        complex_backward_4, complex_backward_4096, complex_backward_64, complex_backward_8,
    };
    use crate::split_radix::complex_forward::{
        complex_forward_128, complex_forward_16, complex_forward_256, complex_forward_32,
        complex_forward_4, complex_forward_4096, complex_forward_64, complex_forward_8,
    };
    use crate::split_radix::{reduce_2p_sqr, Complex, Real, P};

    impl Add for Complex {
        type Output = Self;

        fn add(self, other: Self) -> Self {
            Self {
                re: reduce_2p_sqr(self.re + other.re),
                im: reduce_2p_sqr(self.im + other.im),
            }
        }
    }

    impl Mul for Complex {
        type Output = Self;

        fn mul(self, other: Self) -> Self {
            Self {
                re: reduce_2p_sqr(
                    (self.re % P * (other.re % P)) % P - (self.im % P * (other.im % P)) % P,
                ),
                im: reduce_2p_sqr(
                    (self.im % P * (other.re % P)) % P + (self.re % P * (other.im % P)) % P,
                ),
            }
        }
    }

    impl Mul<Real> for Complex {
        type Output = Self;

        fn mul(self, other: Real) -> Self {
            Self {
                re: reduce_2p_sqr(self.re * other),
                im: reduce_2p_sqr(self.im * other),
            }
        }
    }

    /*
    fn naive_convolve(x: &Vec<Complex>, y: &Vec<Complex>) -> Vec<Complex> {
        let n = x.len();
        assert_eq!(n, y.len());
        let mut z = vec![Complex::ZERO; n];
        for i in 0..n {
            for j in 0..n {
                let mut t = z[(i + j) % n];
                t = t + x[i] * y[j];
                z[(i + j) % n] = t;
            }
        }
        z
    }
    */

    fn naive_convolve(us: &[Complex], vs: &[Complex]) -> Vec<Complex> {
        let n = us.len();
        assert_eq!(n, vs.len());

        let mut conv = Vec::with_capacity(n);
        for i in 0..n {
            let mut t = Complex::ZERO;
            for j in 0..n {
                t = t + us[j] * vs[(n + i - j) % n];
            }
            conv.push(t);
        }
        conv
    }

    fn randcomplex() -> Complex {
        let mut rng = thread_rng();
        let re = rng.gen::<u32>();
        let im = rng.gen::<u32>();
        Complex::new(reduce_2p_sqr(re as i64), reduce_2p_sqr(im as i64))
    }

    fn randvec(n: usize) -> Vec<Complex> {
        repeat_with(randcomplex).take(n).collect::<Vec<_>>()
    }

    #[test]
    fn forward_backward_is_identity4() {
        const N: usize = 4;

        let us = randvec(N);
        let mut vs = us.clone();
        complex_forward_4(&mut vs);

        let mut ws = vs.clone();
        complex_backward_4(&mut ws);

        assert!(us.iter().zip(ws).all(|(&u, w)| w == u * N as i64));
    }

    #[test]
    fn forward_backward_is_identity8() {
        for _ in 0..1000 {
            const N: usize = 8;

            let us = randvec(N);
            let mut vs = us.clone();
            complex_forward_8(&mut vs);

            let mut ws = vs.clone();
            complex_backward_8(&mut ws);

            assert!(us.iter().zip(ws).all(|(&u, w)| w == u * N as i64));
        }
    }

    #[test]
    fn forward_backward_is_identity16() {
        for _ in 0..1000 {
            const N: usize = 16;

            let us = randvec(N);
            let mut vs = us.clone();
            complex_forward_16(&mut vs);

            let mut ws = vs.clone();
            complex_backward_16(&mut ws);

            assert!(us.iter().zip(ws).all(|(&u, w)| w == u * N as i64));
        }
    }

    #[test]
    fn forward_backward_is_identity32() {
        for _ in 0..1000 {
            const N: usize = 32;

            let us = randvec(N);
            let mut vs = us.clone();
            complex_forward_32(&mut vs);

            let mut ws = vs.clone();
            complex_backward_32(&mut ws);

            assert!(us.iter().zip(ws).all(|(&u, w)| w == u * N as i64));
        }
    }

    #[test]
    fn forward_backward_is_identity64() {
        for _ in 0..1000 {
            const N: usize = 64;

            let us = randvec(N);
            let mut vs = us.clone();
            complex_forward_64(&mut vs);

            let mut ws = vs.clone();
            complex_backward_64(&mut ws);

            assert!(us.iter().zip(ws).all(|(&u, w)| w == u * N as i64));
        }
    }

    #[test]
    fn forward_backward_is_identity128() {
        for _ in 0..1000 {
            const N: usize = 128;

            let us = randvec(N);
            let mut vs = us.clone();
            complex_forward_128(&mut vs);

            let mut ws = vs.clone();
            complex_backward_128(&mut ws);

            assert!(us.iter().zip(ws).all(|(&u, w)| w == u * N as i64));
        }
    }

    #[test]
    fn forward_backward_is_identity256() {
        for _ in 0..1000 {
            const N: usize = 256;

            let us = randvec(N);
            let mut vs = us.clone();
            complex_forward_256(&mut vs);

            let mut ws = vs.clone();
            complex_backward_256(&mut ws);

            assert!(us.iter().zip(ws).all(|(&u, w)| w == u * N as i64));
        }
    }

    #[test]
    fn forward_backward_is_identity4096() {
        for _ in 0..100 {
            const N: usize = 4096;

            let us = randvec(N);
            let mut vs = us.clone();
            complex_forward_4096(&mut vs);

            let mut ws = vs.clone();
            complex_backward_4096(&mut ws);

            assert!(us.iter().zip(ws).all(|(&u, w)| w == u * N as i64));
        }
    }

    #[test]
    fn convolution16() {
        const N: usize = 16;

        let us = randvec(N);
        let vs = randvec(N);

        let mut fft_us = us.clone();
        complex_forward_16(&mut fft_us);

        let mut fft_vs = vs.clone();
        complex_forward_16(&mut fft_vs);

        let mut pt_prods = fft_us
            .iter()
            .zip(fft_vs)
            .map(|(&u, v)| u * v)
            .collect::<Vec<_>>();

        complex_backward_16(&mut pt_prods);

        let conv = naive_convolve(&us, &vs);

        assert!(conv.iter().zip(pt_prods).all(|(&c, p)| p == c * N as i64));
    }

    #[test]
    fn convolution32() {
        const N: usize = 32;

        let us = randvec(N);
        let vs = randvec(N);

        let mut fft_us = us.clone();
        complex_forward_32(&mut fft_us);

        let mut fft_vs = vs.clone();
        complex_forward_32(&mut fft_vs);

        let mut pt_prods = fft_us
            .iter()
            .zip(fft_vs)
            .map(|(&u, v)| u * v)
            .collect::<Vec<_>>();

        complex_backward_32(&mut pt_prods);

        let conv = naive_convolve(&us, &vs);

        assert!(conv.iter().zip(pt_prods).all(|(&c, p)| p == c * N as i64));
    }

    #[test]
    fn convolution4096() {
        const N: usize = 4096;

        let us = randvec(N);
        let vs = randvec(N);

        let mut fft_us = us.clone();
        complex_forward_4096(&mut fft_us);

        let mut fft_vs = vs.clone();
        complex_forward_4096(&mut fft_vs);

        let mut pt_prods = fft_us
            .iter()
            .zip(fft_vs)
            .map(|(&u, v)| u * v)
            .collect::<Vec<_>>();

        complex_backward_4096(&mut pt_prods);

        let conv = naive_convolve(&us, &vs);

        assert!(conv.iter().zip(pt_prods).all(|(&c, p)| p == c * N as i64));
    }

    /*
    #[test]
    fn bad_reduc() {
        // In particular -906276279 + 967747991 i is a 32nd root of
        // unity and
        // Re[(-906276279 + 967747991 i)*(1355777538 + 3495721022 i)] = -2^62.

        // These satisfy the ranges as 906276279 < 967747991 <
        // 1073741823 = (P - 1)/2 and 1355777538 < 3495721022 <
        // 4294967294 = 2P.

        let a = Complex::new(-906276279, 967747991); // 32nd root of 1
        let b = Complex::new(1355777538, 3495721022);
        let ab_re = a.re * b.re - a.im * b.im;
        assert_eq!(ab_re, -1 << 62);
        assert_eq!(normalise_real(ab_re), P - 1);
    }
    */
}
