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

use crate::split_radix::{complex_backward::*, complex_forward::*};
use core::fmt;

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
        16 => c16(u),
        32 => c32(u),
        64 => c64(u),
        512 => c512(u),
        1024 => c1024(u),
        2048 => c2048(u),
        4096 => c4096(u),
        _ => panic!("unsupported size"),
    }
}

pub fn backward_fft<const N: usize>(v: &mut [Complex]) {
    match N {
        512 => backward512(v),
        1024 => backward1024(v),
        2048 => backward2048(v),
        4096 => backward4096(v),
        _ => panic!("unsupported size"),
    }
}

const P: Real = (1 << 31) - 1;

impl Complex {
    const ZERO: Self = Self::new(0i64, 0i64);
    //const ONE: Self = Self::new(1i64, 0i64);
    //const J: Self = Self::new(0i64, 1i64);

    pub const fn new(re: Real, im: Real) -> Self {
        Self { re, im }
    }

    pub const fn from_pair((re, im): (Real, Real)) -> Self {
        Self::new(re, im)
    }

    pub fn conj(&self) -> Self {
        Self::new(self.re, normalise_real(-self.im))
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

#[inline]
pub(crate) fn normalise_u32(x: Real) -> Real {
    debug_assert!(x >= 0 && x <= 2 * P);

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

/// Given a Real x in the range |x| <= 2 * P * (P - 1), return
/// the Real x' in the range 0 <= x' <= P such that x' = x (mod P).
#[inline]
pub(crate) fn normalise_real(x: Real) -> Real {
    const ABSMAX_ELT: i64 = 2 * P * (P - 1);
    debug_assert!(x <= ABSMAX_ELT && x >= -ABSMAX_ELT);

    const MASK: u64 = (1 << 31) - 1;
    let y = x as u64; // make unsigned so the right shift doesn't sign extend
    let (lo, hi) = ((y & MASK) as i64, ((y >> 31) & MASK) as i64);

    // Here we prove the correctness of the line "lo + hi + (x >> 62)".
    //
    // FIXME: Note that for one particular value (3 << 62) we "underflow"
    // and pass -1 to `normalise_u32()`. This is bug.
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
    //     This is the signed 2-bit number in 2's complement, so can be
    //     obtained as the arithmetic right shift of the top two bits,
    //     whence (x >> 62) = -2 * b[63] + b[62].
    //  3) Regarding overflow and underflow, consider each possible value
    //     of -2 * b[63] + b[62]:
    //       a.  0 ==> Then the result is lo + hi which is fine.
    //       b.  1 ==> Would be > 2P if lo = hi = P, however that corresponds
    //              to the bit pattern b[63] = 0 and b[i] = 1 for 0 <= i < 63
    //              which is 2^63 - 1 > 2 * P * (P - 1), hence can't happen.
    //       c. -2 ==> Would be < 0 if lo + hi = 0 or 1. These possibilites
    //              correspond to x being:
    //                 0x8000_0000_0000_0000   lo = hi = 0
    //                 0x8000_0000_0000_0001   lo = 1, hi = 0
    //                 0x8000_0000_8000_0000   lo = 0, hi = 1
    //              all of which are < -2 * P * (P - 1), hence can't happen.
    //       d. -1 ==> Will be < 0 if lo = hi = 0. This corresponds to
    //              x = 0xC000_0000_0000_0000 > -2 * P * (P - 1), which
    //              is a valid input (specifically it is = -1 (mod p)).
    //              *** HENCE THIS IS A BUG!! ***
    //
    // x must be signed for the SAR here
    let fold = (lo + hi) + (x >> 62);
    // 0 <= fold <= 2P

    let x_mod_p = normalise_u32(fold);
    // 0 <= x_mod_p <= P

    x_mod_p
}

#[inline]
pub(crate) fn normalise(z: &mut Complex) {
    z.re = normalise_real(z.re);
    z.im = normalise_real(z.im);
}

#[inline]
pub(crate) fn normalise_all(zs: &mut [Complex]) {
    zs.iter_mut().for_each(normalise);
}

impl PartialEq for Complex {
    fn eq(&self, other: &Self) -> bool {
        normalise_real(self.re) == normalise_real(other.re)
            && normalise_real(self.im) == normalise_real(other.im)
    }
}

#[cfg(test)]
mod tests {
    use crate::split_radix::complex_backward::{
        backward128, backward256, backward32, backward4096, backward64, u16, u4, u8,
    };
    use crate::split_radix::complex_forward::{c128, c16, c256, c32, c4, c4096, c64, c8};
    use crate::split_radix::{normalise_real, Complex, Real, P};
    use rand::{thread_rng, Rng};

    use alloc::vec::Vec;
    use core::iter::repeat_with;
    use core::ops::{Add, Mul};

    impl Add for Complex {
        type Output = Self;

        fn add(self, other: Self) -> Self {
            Self {
                re: normalise_real(self.re + other.re),
                im: normalise_real(self.im + other.im),
            }
        }
    }

    impl Mul for Complex {
        type Output = Self;

        fn mul(self, other: Self) -> Self {
            Self {
                re: normalise_real(
                    (self.re % P * (other.re % P)) % P - (self.im % P * (other.im % P)) % P,
                ),
                im: normalise_real(
                    (self.im % P * (other.re % P)) % P + (self.re % P * (other.im % P)) % P,
                ),
            }
        }
    }

    impl Mul<Real> for Complex {
        type Output = Self;

        fn mul(self, other: Real) -> Self {
            Self {
                re: normalise_real(self.re * other),
                im: normalise_real(self.im * other),
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
        Complex::new(normalise_real(re as i64), normalise_real(im as i64))
    }

    fn randvec(n: usize) -> Vec<Complex> {
        repeat_with(randcomplex).take(n).collect::<Vec<_>>()
    }

    #[test]
    fn forward_backward_is_identity4() {
        const N: usize = 4;

        let us = randvec(N);
        let mut vs = us.clone();
        c4(&mut vs);

        let mut ws = vs.clone();
        u4(&mut ws);

        assert!(us.iter().zip(ws).all(|(&u, w)| w == u * N as i64));
    }

    #[test]
    fn forward_backward_is_identity8() {
        for _ in 0..1000 {
            const N: usize = 8;

            let us = randvec(N);
            let mut vs = us.clone();
            c8(&mut vs);

            let mut ws = vs.clone();
            u8(&mut ws);

            assert!(us.iter().zip(ws).all(|(&u, w)| w == u * N as i64));
        }
    }

    #[test]
    fn forward_backward_is_identity16() {
        for _ in 0..1000 {
            const N: usize = 16;

            let us = randvec(N);
            let mut vs = us.clone();
            c16(&mut vs);

            let mut ws = vs.clone();
            u16(&mut ws);

            assert!(us.iter().zip(ws).all(|(&u, w)| w == u * N as i64));
        }
    }

    #[test]
    fn forward_backward_is_identity32() {
        for _ in 0..1000 {
            const N: usize = 32;

            let us = randvec(N);
            let mut vs = us.clone();
            c32(&mut vs);

            let mut ws = vs.clone();
            backward32(&mut ws);

            assert!(us.iter().zip(ws).all(|(&u, w)| w == u * N as i64));
        }
    }

    #[test]
    fn forward_backward_is_identity64() {
        for _ in 0..1000 {
            const N: usize = 64;

            let us = randvec(N);
            let mut vs = us.clone();
            c64(&mut vs);

            let mut ws = vs.clone();
            backward64(&mut ws);

            assert!(us.iter().zip(ws).all(|(&u, w)| w == u * N as i64));
        }
    }

    #[test]
    fn forward_backward_is_identity128() {
        for _ in 0..1000 {
            const N: usize = 128;

            let us = randvec(N);
            let mut vs = us.clone();
            c128(&mut vs);

            let mut ws = vs.clone();
            backward128(&mut ws);

            assert!(us.iter().zip(ws).all(|(&u, w)| w == u * N as i64));
        }
    }

    #[test]
    fn forward_backward_is_identity256() {
        for _ in 0..1000 {
            const N: usize = 256;

            let us = randvec(N);
            let mut vs = us.clone();
            c256(&mut vs);

            let mut ws = vs.clone();
            backward256(&mut ws);

            assert!(us.iter().zip(ws).all(|(&u, w)| w == u * N as i64));
        }
    }

    #[test]
    fn forward_backward_is_identity4096() {
        for _ in 0..100 {
            const N: usize = 4096;

            let us = randvec(N);
            let mut vs = us.clone();
            c4096(&mut vs);

            let mut ws = vs.clone();
            backward4096(&mut ws);

            assert!(us.iter().zip(ws).all(|(&u, w)| w == u * N as i64));
        }
    }

    #[test]
    fn convolution16() {
        const N: usize = 16;

        let us = randvec(N);
        let vs = randvec(N);

        let mut fft_us = us.clone();
        c16(&mut fft_us);

        let mut fft_vs = vs.clone();
        c16(&mut fft_vs);

        let mut pt_prods = fft_us
            .iter()
            .zip(fft_vs)
            .map(|(&u, v)| u * v)
            .collect::<Vec<_>>();

        u16(&mut pt_prods);

        let conv = naive_convolve(&us, &vs);

        assert!(conv.iter().zip(pt_prods).all(|(&c, p)| p == c * N as i64));
    }

    #[test]
    fn convolution32() {
        const N: usize = 32;

        let us = randvec(N);
        let vs = randvec(N);

        let mut fft_us = us.clone();
        c32(&mut fft_us);

        let mut fft_vs = vs.clone();
        c32(&mut fft_vs);

        let mut pt_prods = fft_us
            .iter()
            .zip(fft_vs)
            .map(|(&u, v)| u * v)
            .collect::<Vec<_>>();

        backward32(&mut pt_prods);

        let conv = naive_convolve(&us, &vs);

        assert!(conv.iter().zip(pt_prods).all(|(&c, p)| p == c * N as i64));
    }

    #[test]
    fn convolution4096() {
        const N: usize = 4096;

        let us = randvec(N);
        let vs = randvec(N);

        let mut fft_us = us.clone();
        c4096(&mut fft_us);

        let mut fft_vs = vs.clone();
        c4096(&mut fft_vs);

        let mut pt_prods = fft_us
            .iter()
            .zip(fft_vs)
            .map(|(&u, v)| u * v)
            .collect::<Vec<_>>();

        backward4096(&mut pt_prods);

        let conv = naive_convolve(&us, &vs);

        assert!(conv.iter().zip(pt_prods).all(|(&c, p)| p == c * N as i64));
    }
}
