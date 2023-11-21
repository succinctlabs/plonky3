use super::roots::{D1024, D128, D16, D2048, D256, D32, D4096, D512, D64, D8192};
use super::{reduce_2p, reduce_complex, split_at_mut_unchecked, Complex, Real, P};

/// SQRTHALF * SQRTHALF = 1/2 (mod P)
const SQRTHALF: Real = D16[1].re; // == 1 << 15

// a0 = (a0.re + a2.re, a0.im + a2.im)        = a0 + a2
// a1 = (a1.re + a3.re, a1.im + a3.im)        = a1 + a3
// a2 = (
//        wre * [(a0.re - a2.re) - (a1.im - a3.im)]
//        - wim * [(a0.im - a2.im) + (a1.re - a3.re)]
//       ,
//        wre * [(a0.im - a2.im) + (a1.re - a3.re)]
//        + wim * [(a0.re - a2.re) - (a1.im - a3.im)]
//      )
//
// a3 = (
//        wre * [(a0.re - a2.re) + (a1.im - a3.im)]
//        + wim * [(a0.im - a2.im) - (a1.re - a3.re)]
//       ,
//        wre * [(a0.im - a2.im) - (a1.re - a3.re)]
//        - wim * [(a0.re - a2.re) + (a1.im - a3.im)]
//      )
//
// Input/output are Complex a0, a1, a2, a3. Assume inputs are
// "essentially reduced" modulo P (i.e. 0 <= ai.{re,im} <= P).
//
// Root (wre, wim) is in "balanced representation", i.e.
// -(P-1)/2 <= w.{re,im} <= (P-1)/2.
//
// Then
//
//   0 <= a0.{re,im} <= 2P
//   0 <= a1.{re,im} <= 2P
//
// NB: These inequalities can't all attain a max (resp. min) simultaneously.
//
// |a[23].{re,im}| <= 2 * 2P * (P - 1)/2 = 2P^2 - 2P < 2^63
//
// For reference, if the roots were in reduced form, then the bounds would be:
//
// |a[23].{re,im}| <= 2 * ((P - 1) * 2P) = 4P^2 - 4P > 2^64
//
#[inline]
#[rustfmt::skip]
fn transform(
    a0: &mut Complex,
    a1: &mut Complex,
    a2: &mut Complex,
    a3: &mut Complex,
    wre: Real,
    wim: Real,
) {
    debug_assert!(0 <= a0.re, "{}", a0.re);
    debug_assert!(a0.re <= P, "{}", a0.re);
    debug_assert!(0 <= a0.im, "{}", a0.im);
    debug_assert!(a0.im <= P, "{}", a0.im);

    debug_assert!(0 <= a1.re);
    debug_assert!(a1.re <= P);
    debug_assert!(0 <= a1.im);
    debug_assert!(a1.im <= P);

    debug_assert!(0 <= a2.re);
    debug_assert!(a2.re <= P);
    debug_assert!(0 <= a2.im);
    debug_assert!(a2.im <= P);

    debug_assert!(0 <= a3.re);
    debug_assert!(a3.re <= P);
    debug_assert!(0 <= a3.im);
    debug_assert!(a3.im <= P);

    let mut t6 = a2.re;       // t6 = a2.re
    let mut t1 = a0.re - t6;  // t1 = a0.re - a2.re
    t6 += a0.re;              // t6 = a0.re + a2.re
    a0.re = reduce_2p(t6);
    let mut t3 = a3.im;       // t3 = a3.im
    let mut t4 = a1.im - t3;  // t4 = a1.im - a3.im
    let mut t8 = t1 - t4;     // t8 = (a0.re - a2.re) - (a1.im - a3.im)
    t1 += t4;                 // t1 = (a0.re - a2.re) + (a1.im - a3.im)
    t3 += a1.im;              // t3 = a1.im + a3.im
    a1.im = reduce_2p(t3);
    let mut t5 = wre;         // t5 = wre
    let mut t7 = t8 * t5;     // t7 = wre * [(a0.re - a2.re) - (a1.im - a3.im)]
    t4 = t1 * t5;             // t4 = wre * [(a0.re - a2.re) + (a1.im - a3.im)]
    t8 *= wim;                // t8 = wim * [(a0.re - a2.re) - (a1.im - a3.im)]
    let mut t2 = a3.re;       // t2 = a3.re
    t3 = a1.re - t2;          // t3 = a1.re - a3.re
    t2 += a1.re;              // t2 = a1.re + a3.re
    a1.re = reduce_2p(t2);

    t1 *= wim;                // t1 = wim * [(a0.re - a2.re) + (a1.im - a3.im)]
    t6 = a2.im;               // t6 = a2.im
    t2 = a0.im - t6;          // t2 = a0.im - a2.im
    t6 += a0.im;              // t6 = a0.im + a2.im
    a0.im = reduce_2p(t6);

    t6 = t2 + t3;             // t6 = (a0.im - a2.im) + (a1.re - a3.re)
    t2 -= t3;                 // t2 = (a0.im - a2.im) - (a1.re - a3.re)
    t3 = t6 * wim;            // t3 = wim * [(a0.im - a2.im) + (a1.re - a3.re)]
    t7 -= t3;                 // t7 = wre * [(a0.re - a2.re) - (a1.im - a3.im)]
                              //     - wim * [(a0.im - a2.im) + (a1.re - a3.re)]
    a2.re = t7;
    t6 *= t5;                 // t6 = wre * [(a0.im - a2.im) + (a1.re - a3.re)]
    t6 += t8;                 // t6 = wre * [(a0.im - a2.im) + (a1.re - a3.re)]
                              //     + wim * [(a0.re - a2.re) - (a1.im - a3.im)]
    a2.im = t6;

    reduce_complex(a2);

    t5 *= t2;                 // t5 = wre * [(a0.im - a2.im) - (a1.re - a3.re)]
    t5 -= t1;                 // t5 = wre * [(a0.im - a2.im) - (a1.re - a3.re)]
                              //     - wim * [(a0.re - a2.re) + (a1.im - a3.im)]
    a3.im = t5;
    t2 *= wim;                // t2 = wim * [(a0.im - a2.im) - (a1.re - a3.re)]
    t4 += t2;                 // t4 = wre * [(a0.re - a2.re) + (a1.im - a3.im)]
                              //     + wim * [(a0.im - a2.im) - (a1.re - a3.re)]
    a3.re = t4;

    reduce_complex(a3);
}

// Appears to be related to applying the radix-8 butterfly
//
// a0 = (a0.re + a2.re, a0.im + a2.im)        = a0 + a2
// a1 = (a1.re + a3.re, a1.im + a3.im)        = a1 + a3
// a2 = ([(a0.re - a2.re) - (a1.im - a3.im)
//          - ((a0.im - a2.im) + (a1.re - a3.re))] * sqrthalf,
//       [(a0.re - a2.re) - (a1.im - a3.im)
//          + ((a0.im - a2.im) + (a1.re - a3.re))] * sqrthalf)
// a3 = ([(a0.re - a2.re) + (a1.im - a3.im)
//          + ((a0.im - a2.im) - (a1.re - a3.re))] * sqrthalf,
//       [-((a0.re - a2.re) + (a1.im - a3.im))
//          + ((a0.im - a2.im) - (a1.re - a3.re))] * sqrthalf)
//
// Input/output are Complex a0, a1, a2, a3. Assume inputs are
// "essentially reduced" modulo P (i.e. 0 <= ai.{re,im} <= P).
//
// NB: sqrthalf = 2^15.
//
// Then
//
//   0 <= a0.{re,im} <= 2P
//   0 <= a1.{re,im} <= 2P
//
// NB: These inequalities can't all attain a max (resp. min) simultaneously.
//
// -2^17 P <= a[23].{re,im} <= 2^17 P
//
#[inline]
fn transformhalf(a0: &mut Complex, a1: &mut Complex, a2: &mut Complex, a3: &mut Complex) {
    let mut t1 = a2.re;
    let mut t5 = a0.re - t1;
    t1 += a0.re;
    a0.re = reduce_2p(t1);
    let mut t4 = a3.im;
    let mut t8 = a1.im - t4;
    t1 = t5 - t8;
    t5 += t8;
    t4 += a1.im;
    a1.im = reduce_2p(t4);
    let mut t3 = a3.re;
    let mut t7 = a1.re - t3;
    t3 += a1.re;
    a1.re = reduce_2p(t3);

    t8 = a2.im;
    let mut t6 = a0.im - t8;
    let mut t2 = t6 + t7;
    t6 -= t7;
    t8 += a0.im;
    a0.im = reduce_2p(t8);

    t4 = t6 + t5;
    t3 = SQRTHALF;
    t4 *= t3;
    a3.re = t4;
    t6 -= t5;
    t6 *= t3;
    a3.im = t6;

    // TODO: These reductions mightn't be necessary, as removing them
    // doesn't trigger overflow in the tests.
    reduce_complex(a3);

    t7 = t1 - t2;
    t7 *= t3;
    a2.re = t7;
    t2 += t1;
    t2 *= t3;
    a2.im = t2;

    reduce_complex(a2);
}

// a0 = (a0.re + a2.re, a0.im + a2.im)        = a0 + a2
// a1 = (a1.re + a3.re, a1.im + a3.im)        = a1 + a3
// a2 = ((a0.re - a2.re) - (a1.im - a3.im),
//       (a0.im - a2.im) + (a1.re - a3.re))   = (a0 - a2) + i*(a1 - a3)
// a3 = ((a0.re - a2.re) + (a1.im - a3.im),
//       (a0.im - a2.im) - (a1.re - a3.re))   = (a0 - a2) - i*(a1 - a3)
//
// Input/output are Complex a0, a1, a2, a3. Assume inputs are
// "essentially reduced" modulo P (i.e. 0 <= ai.{re,im} <= P).
//
// Then
//
//   0 <= a0.{re,im} <= 2P
//   0 <= a1.{re,im} <= 2P
// -2P <= a2.re, a3.re <= 2P  // NB: These two lines can't attain a max
// -2P <= a2.im, a3.im <= 2P  // (resp. min) simultaneously.
//
#[inline]
fn transformzero(a0: &mut Complex, a1: &mut Complex, a2: &mut Complex, a3: &mut Complex) {
    let mut t5 = a2.re;
    let mut t1 = a0.re - t5;
    t5 += a0.re;
    a0.re = reduce_2p(t5);
    let mut t8 = a3.im;
    let t4 = a1.im - t8;
    let mut t7 = a3.re;
    let mut t6 = t1 - t4;
    a2.re = t6;
    t1 += t4;
    a3.re = t1;
    t8 += a1.im;
    a1.im = reduce_2p(t8);
    let t3 = a1.re - t7;
    t7 += a1.re;
    a1.re = reduce_2p(t7);

    t6 = a2.im;
    let mut t2 = a0.im - t6;
    t7 = t2 + t3;
    a2.im = t7;

    // TODO: Find a cheaper reduction since we know
    // -2P <= a2.{re,im} <= 2P
    reduce_complex(a2);

    t2 -= t3;
    a3.im = t2;

    // TODO: Find a cheaper reduction since we know
    // -2P <= a2.{re,im} <= 2P
    reduce_complex(a3);

    t6 += a0.im;
    a0.im = reduce_2p(t6);
}

// Conjugate pair radix-4 butterfly with twiddle = w = -i (note that
// w^-1 = \bar{w} (the conjugate) when w has unit length):
//
//   s  <-- w a2 + w^-1 a3
//   t  <-- w a2 - w^-1 a3
//   b0 <-- a0 + s    = a0 - i a2 + i a3
//   b1 <-- a1 - i t  = a1 - a2 - a3
//   b2 <-- a2 - s    = a0 + i a2 - i a3
//   b3 <-- a3 + i t  = a1 + a2 + a3
//
// So the transformation matrix is:
//
//   [ 1  0 -i  i ]
//   [ 0  1 -1 -1 ]
//   [ 1  0  i -i ]
//   [ 0  1  1  1 ]
//
// compared with the usual matrix:
//
//   [ 1  1  1  1 ]
//   [ 1 -i -1  i ]
//   [ 1 -1  1 -1 ]
//   [ 1  i -1 -i ]
//

// Expansion of splitfft for n = 4
//
// x[0..4]
//
// u[0..2] = fft(x[0], x[2])
// u0 = x0 + x2
// u1 = x0 - x2
// z  = fft(x[1]) = x1
// z' = fft(x[3]) = x3
//
// w = 1:
// y0 = u0 +   (z + z')  =  x0 + x2 +   x1 +   x3  =  x0 +   x1 + x2 +   x3
// y1 = u1 - i (z - z')  =  x0 - x2 - i x1 + i x3  =  x0 - i x1 - x2 + i x3
// y2 = u0 -   (z + z')  =  x0 + x2 -   x1 -   x3  =  x0 -   x1 + x2 -   x3
// y3 = u1 + i (z - z')  =  x0 - x2 + i x1 - i x3  =  x0 + i x1 - x2 - i x3

// This is the usual matrix (above)!
//
// [ 1  1  1  1 ]
// [ 1 -i -1  i ]
// [ 1 -1  1 -1 ]
// [ 1  i -1 -i ]

// y0 = u0 +   (w z + w^-1 z') = x0 + x2 +  w x1 +  w^-1 x3
// y2 = u0 -   (w z + w^-1 z') = x0 + x2 -  w x1 -  w^-1 x3
// y1 = u1 - i (w z - w^-1 z') = x0 - x2 - iw x1 + iw^-1 x3
// y3 = u1 + i (w z - w^-1 z') = x0 - x2 + iw x1 - iw^-1 x3
//

#[inline(always)]
#[rustfmt::skip]
pub(crate) fn complex_forward_4(a: &mut [Complex]) {
    assert_eq!(a.len(), 4);

    // TODO: We know that the outputs are in the range -3P <= out <= 4P,
    // so we can possibly use a faster reduction algo.
    
    // This is the conjugate transpose of the matrix for u4.
    //
    // [ 1  1  1  1 ]
    // [ 1 -1  1 -1 ]
    // [ 1  i -1 -i ]
    // [ 1 -i -1  i ]
    //
    // b0 = (a0.re + a1.re + a2.re + a3.re) + i (a0.im + a1.im + a2.im + a3.im)
    // b1 = (a0.re - a1.re + a2.re - a3.re) + i (a0.im - a1.im + a2.im - a3.im)
    // b2 = (a0.re - a1.im - a2.re + a3.im) + i (a0.im + a1.re - a2.im - a3.re)
    // b3 = (a0.re + a1.im - a2.re - a3.im) + i (a0.im - a1.re - a2.im + a3.re)

    /*
    let b0 = Complex::new(a[0].re + a[1].re + a[2].re + a[3].re, a[0].im + a[1].im + a[2].im + a[3].im);
    let b1 = Complex::new(a[0].re - a[1].re + a[2].re - a[3].re, a[0].im - a[1].im + a[2].im - a[3].im);
    let b2 = Complex::new(a[0].re - a[1].im - a[2].re + a[3].im, a[0].im + a[1].re - a[2].im - a[3].re);
    let b3 = Complex::new(a[0].re + a[1].im - a[2].re - a[3].im, a[0].im - a[1].re - a[2].im + a[3].re);

    a[0] = b0;
    a[1] = b1;
    a[2] = b2;
    a[3] = b3;
    */
    
    let mut t5 = a[2].re;       // a2.re
    let mut t1 = a[0].re - t5;  // a0.re - a2.re
    let mut t7 = a[3].re;       // a3.re
    t5 += a[0].re;              // a0.re + a2.re
    let t3 = a[1].re - t7;      // a1.re - a3.re
    t7 += a[1].re;              // a1.re + a3.re
    let t8 = t5 + t7;           // a0.re + a1.re + a2.re + a3.re
    a[0].re = t8;
    t5 -= t7;                   // a0.re - a1.re + a2.re - a3.re
    a[1].re = t5;
    let mut t6 = a[2].im;       // a2.im
    let mut t2 = a[0].im - t6;  // a0.im - a2.im
    t6 += a[0].im;              // a0.im + a2.im
    t5 = a[3].im;               // a3.im
    a[2].im = t2 + t3;          // a0.im + a1.re - a2.im - a3.re
    t2 -= t3;                   // a0.im - a1.re - a2.im + a3.re
    a[3].im = t2;
    let t4 = a[1].im - t5;      // a1.im - a3.im
    a[3].re = t1 + t4;          // a0.re + a1.im - a2.re - a3.im
    t1 -= t4;                   // a0.re - a1.im - a2.re + a3.im
    a[2].re = t1;
    t5 += a[1].im;              // a1.im + a3.im
    a[0].im = t6 + t5;          // a0.im + a1.im + a2.im + a3.im
    t6 -= t5;                   // a0.im - a1.im + a2.im - a3.im
    a[1].im = t6;

    // TODO: Not reducing here produces the right answers, but the
    // result still needs to be reduced at the end.
    //
    // reduce_4p(&mut a[0]);
    // reduce_4p(&mut a[1]);
    // reduce_4p(&mut a[2]);
    // reduce_4p(&mut a[3]);
}

#[inline(always)]
pub(crate) fn complex_forward_8(a: &mut [Complex]) {
    assert_eq!(a.len(), 8);

    let mut t7 = a[4].im;
    let mut t4 = a[0].im - t7;
    t7 += a[0].im;
    a[0].im = t7;

    let mut t8 = a[6].re;
    let mut t5 = a[2].re - t8;
    t8 += a[2].re;
    a[2].re = t8;

    t7 = a[6].im;
    a[6].im = t4 - t5;
    t4 += t5;
    a[4].im = t4;

    let mut t6 = a[2].im - t7;
    t7 += a[2].im;
    a[2].im = t7;

    t8 = a[4].re;
    let mut t3 = a[0].re - t8;
    t8 += a[0].re;
    a[0].re = t8;

    a[4].re = t3 - t6;
    t3 += t6;
    a[6].re = t3;

    t7 = a[5].re;
    t3 = a[1].re - t7;
    t7 += a[1].re;
    a[1].re = t7;

    t8 = a[7].im;
    t6 = a[3].im - t8;
    t8 += a[3].im;
    a[3].im = t8;
    let mut t1 = t3 - t6;
    t3 += t6;

    t7 = a[5].im;
    t4 = a[1].im - t7;
    t7 += a[1].im;
    a[1].im = t7;

    t8 = a[7].re;
    t5 = a[3].re - t8;
    t8 += a[3].re;
    a[3].re = t8;

    let mut t2 = t4 - t5;
    t4 += t5;

    t6 = t1 - t4;
    t8 = SQRTHALF;
    t6 *= t8;
    a[5].re = a[4].re - t6;
    t1 += t4;
    t1 *= t8;
    a[5].im = a[4].im - t1;
    t6 += a[4].re;
    a[4].re = t6;
    t1 += a[4].im;
    a[4].im = t1;

    t5 = t2 - t3;
    t5 *= t8;
    a[7].im = a[6].im - t5;
    t2 += t3;
    t2 *= t8;
    a[7].re = a[6].re - t2;
    t2 += a[6].re;
    a[6].re = t2;
    t5 += a[6].im;
    a[6].im = t5;

    complex_forward_4(&mut a[..4]);

    // TODO: SQRTHALF is known to be smaller than a generic element,
    // so we could specialise to a simpler reduction than normalise_real.
    reduce_complex(&mut a[4]);
    reduce_complex(&mut a[5]);
    reduce_complex(&mut a[6]);
    reduce_complex(&mut a[7]);
}

pub fn complex_forward_8_array(a: &mut [Complex; 8]) {
    complex_forward_8(a);
}

#[inline]
pub(crate) fn complex_forward_16(a: &mut [Complex]) {
    assert_eq!(a.len(), 16);

    let (a0, a1) = unsafe { split_at_mut_unchecked(a, 4) };
    let (a1, a2) = unsafe { split_at_mut_unchecked(a1, 4) };
    let (a2, a3) = unsafe { split_at_mut_unchecked(a2, 4) };

    transformzero(&mut a0[0], &mut a1[0], &mut a2[0], &mut a3[0]);
    transformhalf(&mut a0[2], &mut a1[2], &mut a2[2], &mut a3[2]);
    transform(
        &mut a0[1], &mut a1[1], &mut a2[1], &mut a3[1], D16[0].re, D16[0].im,
    );
    transform(
        &mut a0[3], &mut a1[3], &mut a2[3], &mut a3[3], D16[0].im, D16[0].re,
    );

    complex_forward_4(&mut a[8..12]);
    complex_forward_4(&mut a[12..16]);
    complex_forward_8(&mut a[..8]);
}

pub fn complex_forward_16_array(a: &mut [Complex; 16]) {
    complex_forward_16(a);
}

#[inline]
pub(crate) fn complex_forward_32(a: &mut [Complex]) {
    assert_eq!(a.len(), 32);

    let (a0, a1) = unsafe { split_at_mut_unchecked(a, 8) };
    let (a1, a2) = unsafe { split_at_mut_unchecked(a1, 8) };
    let (a2, a3) = unsafe { split_at_mut_unchecked(a2, 8) };

    transformzero(&mut a0[0], &mut a1[0], &mut a2[0], &mut a3[0]);

    transform(
        &mut a0[1], &mut a1[1], &mut a2[1], &mut a3[1], D32[0].re, D32[0].im,
    );

    transform(
        &mut a0[2], &mut a1[2], &mut a2[2], &mut a3[2], D32[1].re, D32[1].im,
    );

    transform(
        &mut a0[3], &mut a1[3], &mut a2[3], &mut a3[3], D32[2].re, D32[2].im,
    );

    transformhalf(&mut a0[4], &mut a1[4], &mut a2[4], &mut a3[4]);

    transform(
        &mut a0[5], &mut a1[5], &mut a2[5], &mut a3[5], D32[2].im, D32[2].re,
    );
    transform(
        &mut a0[6], &mut a1[6], &mut a2[6], &mut a3[6], D32[1].im, D32[1].re,
    );
    transform(
        &mut a0[7], &mut a1[7], &mut a2[7], &mut a3[7], D32[0].im, D32[0].re,
    );

    complex_forward_8(&mut a[16..24]);
    complex_forward_8(&mut a[24..32]);
    complex_forward_16(&mut a[..16]);
}

pub fn complex_forward_32_array(a: &mut [Complex; 32]) {
    complex_forward_32(a);
}

// a[0...8n-1], w[0...2n-2]; n >= 2
//
// TODO: Original comment is as above, but note that w should have
// length 2n-1, as is obvious from the original code, which addresses
// an odd number of elements of w.
fn complex_forward_pass(a: &mut [Complex], w: &[Complex]) {
    debug_assert_eq!(a.len() % 8, 0);

    let n = a.len() / 8;

    debug_assert!(n >= 2);
    debug_assert_eq!(w.len(), 2 * n - 1);

    // Split a into four chunks of size 2*n.

    let (a0, a1) = unsafe { split_at_mut_unchecked(a, 2 * n) };
    let (a1, a2) = unsafe { split_at_mut_unchecked(a1, 2 * n) };
    let (a2, a3) = unsafe { split_at_mut_unchecked(a2, 2 * n) };

    transformzero(&mut a0[0], &mut a1[0], &mut a2[0], &mut a3[0]);

    // NB: The original version pulled the first iteration out of the
    // loop and unrolled the loop two iterations. When I did that all
    // the larger (>32) FFTs slowed down, most by about 25-35%!
    for i in 1..2 * n {
        transform(
            &mut a0[i],
            &mut a1[i],
            &mut a2[i],
            &mut a3[i],
            w[i - 1].re,
            w[i - 1].im,
        );
    }
}

fn complex_forward_pass_big(a: &mut [Complex], w: &[Complex]) {
    debug_assert_eq!(a.len() % 8, 0);

    let n = a.len() / 8;

    debug_assert!(n >= 2);
    debug_assert_eq!(w.len(), 2 * n - 1);

    // Split a into four chunks of size 2*n.

    let (a0, a1) = unsafe { split_at_mut_unchecked(a, 2 * n) };
    let (a1, a2) = unsafe { split_at_mut_unchecked(a1, 2 * n) };
    let (a2, a3) = unsafe { split_at_mut_unchecked(a2, 2 * n) };

    transformzero(&mut a0[0], &mut a1[0], &mut a2[0], &mut a3[0]);
    for i in 1..n {
        transform(
            &mut a0[i],
            &mut a1[i],
            &mut a2[i],
            &mut a3[i],
            w[i - 1].re,
            w[i - 1].im,
        );
    }

    transformhalf(&mut a0[n], &mut a1[n], &mut a2[n], &mut a3[n]);

    for i in (n + 1)..(2 * n) {
        transform(
            &mut a0[i],
            &mut a1[i],
            &mut a2[i],
            &mut a3[i],
            w[2 * n - i - 1].im,
            w[2 * n - i - 1].re,
        );
    }
}

// Unrolling this as for c32 had some minor gains for N <= 128, but
// minor slow-downs for bigger sizes.
#[inline]
pub(crate) fn complex_forward_64(a: &mut [Complex]) {
    assert_eq!(a.len(), 64);

    complex_forward_pass(a, &D64);
    complex_forward_16(&mut a[32..48]);
    complex_forward_16(&mut a[48..64]);
    complex_forward_32(&mut a[..32]);
}

pub fn complex_forward_64_array(a: &mut [Complex; 64]) {
    complex_forward_64(a);
}

pub(crate) fn complex_forward_128(a: &mut [Complex]) {
    assert_eq!(a.len(), 128);

    complex_forward_pass(a, &D128);
    complex_forward_32(&mut a[64..96]);
    complex_forward_32(&mut a[96..128]);
    complex_forward_64(&mut a[..64]);
}

pub(crate) fn complex_forward_256(a: &mut [Complex]) {
    assert_eq!(a.len(), 256);

    complex_forward_pass(a, &D256);
    complex_forward_64(&mut a[128..192]);
    complex_forward_64(&mut a[192..256]);
    complex_forward_128(&mut a[..128]);
}

pub(crate) fn complex_forward_512(a: &mut [Complex]) {
    assert_eq!(a.len(), 512);

    complex_forward_pass(a, &D512);
    complex_forward_128(&mut a[384..512]);
    complex_forward_128(&mut a[256..384]);
    complex_forward_256(&mut a[..256]);
}

pub(crate) fn complex_forward_1024(a: &mut [Complex]) {
    assert_eq!(a.len(), 1024);

    complex_forward_pass_big(a, &D1024);
    complex_forward_256(&mut a[768..1024]);
    complex_forward_256(&mut a[512..768]);
    complex_forward_512(&mut a[..512]);
}

pub(crate) fn complex_forward_2048(a: &mut [Complex]) {
    assert_eq!(a.len(), 2048);

    complex_forward_pass_big(a, &D2048);
    complex_forward_512(&mut a[1536..2048]);
    complex_forward_512(&mut a[1024..1536]);
    complex_forward_1024(&mut a[..1024]);
}

pub(crate) fn complex_forward_4096(a: &mut [Complex]) {
    assert_eq!(a.len(), 4096);

    complex_forward_pass_big(a, &D4096);
    complex_forward_1024(&mut a[3072..4096]);
    complex_forward_1024(&mut a[2048..3072]);
    complex_forward_2048(&mut a[..2048]);
}

pub(crate) fn complex_forward_8192(a: &mut [Complex]) {
    assert_eq!(a.len(), 8192);

    complex_forward_pass_big(a, &D8192);
    complex_forward_2048(&mut a[6144..8192]);
    complex_forward_2048(&mut a[4096..6144]);
    complex_forward_4096(&mut a[..4096]);
}
