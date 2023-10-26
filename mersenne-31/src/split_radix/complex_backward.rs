use super::roots::{D1024, D128, D16, D2048, D256, D32, D4096, D512, D64};
use super::{normalise, normalise_all, normalise_i64, Complex, Real};

/// SQRTHALF * SQRTHALF = 1/2 (mod P)
const SQRTHALF: Real = D16[1].re; // == 1 << 15

fn untransform(
    a0: &mut Complex,
    a1: &mut Complex,
    a2: &mut Complex,
    a3: &mut Complex,
    wre: Real,
    wim: Real,
) {
    let mut t6 = wre;
    let mut t1 = a2.re;
    t1 *= t6;
    let mut t8 = wim;
    let mut t3 = a2.im;
    t3 *= t8;
    let mut t2 = a2.im;
    let mut t4 = a2.re;
    let mut t5 = a3.re;
    t5 *= t6;
    t5 = normalise_i64(t5);
    let mut t7 = a3.im;
    t1 += t3;
    t1 = normalise_i64(t1);
    t7 *= t8;
    t5 -= t7;
    t3 = t5 + t1;
    t5 -= t1;
    t2 *= t6;
    t6 *= a3.im;
    t6 = normalise_i64(t6);
    t4 *= t8;
    t2 -= t4;
    t8 *= a3.re;
    t6 += t8;
    t1 = a0.re - t3;
    t3 += a0.re;
    a0.re = t3;
    t7 = a1.im - t5;
    t5 += a1.im;
    a1.im = t5;
    t4 = t2 - t6;
    t6 += t2;
    t8 = a1.re - t4;
    t4 += a1.re;
    a1.re = t4;
    t2 = a0.im - t6;
    t6 += a0.im;
    a0.im = t6;
    a2.re = t1;
    a3.im = t7;
    a3.re = t8;
    a2.im = t2;

    normalise(a0);
    normalise(a1);
    normalise(a2);
    normalise(a3);
}

fn untransformhalf(a0: &mut Complex, a1: &mut Complex, a2: &mut Complex, a3: &mut Complex) {
    let mut t6 = SQRTHALF;
    let mut t1 = a2.re;
    let mut t2 = a2.im - t1;
    t2 *= t6;
    t1 += a2.im;
    t1 *= t6;
    let mut t4 = a3.im;
    let mut t3 = a3.re - t4;
    t3 *= t6;
    t4 += a3.re;
    t4 *= t6;
    let mut t8 = t3 - t1;
    let mut t7 = t2 - t4;
    t1 += t3;
    t2 += t4;
    t4 = a1.im - t8;
    a3.im = t4;
    t8 += a1.im;
    a1.im = t8;
    t3 = a1.re - t7;
    a3.re = t3;
    t7 += a1.re;
    a1.re = t7;
    let t5 = a0.re - t1;
    a2.re = t5;
    t1 += a0.re;
    a0.re = t1;
    t6 = a0.im - t2;
    a2.im = t6;
    t2 += a0.im;
    a0.im = t2;

    normalise(a0);
    normalise(a1);
    normalise(a2);
    normalise(a3);
}

fn untransformzero(a0: &mut Complex, a1: &mut Complex, a2: &mut Complex, a3: &mut Complex) {
    let mut t2 = a3.im;
    let mut t3 = a2.im - t2;
    t2 += a2.im;
    let mut t1 = a2.re;
    let mut t4 = a3.re - t1;
    t1 += a3.re;
    let t5 = a0.re - t1;
    a2.re = t5;
    let t6 = a0.im - t2;
    a2.im = t6;
    let t7 = a1.re - t3;
    a3.re = t7;
    let t8 = a1.im - t4;
    a3.im = t8;
    t1 += a0.re;
    a0.re = t1;
    t2 += a0.im;
    a0.im = t2;
    t3 += a1.re;
    a1.re = t3;
    t4 += a1.im;
    a1.im = t4;

    normalise(a0);
    normalise(a1);
    normalise(a2);
    normalise(a3);
}

#[inline]
pub(crate) fn u4(a: &mut [Complex]) {
    // This is the usual matrix with the middle columns reversed
    // because of the use of w^-1 in place of w^3k.
    //
    // It is the conjugate transpose of the matrix for c4, and we also
    // have  D = 4 * C^-1  (where C is the c4 matrix and D is this matrix)
    //
    // [ 1  1  1  1 ]
    // [ 1 -1 -i  i ]
    // [ 1  1 -1 -1 ]
    // [ 1 -1  i -i ]
    //
    // b0 = (a0.re + a1.re + a2.re + a3.re) + i (a0.im + a1.im + a2.im + a3.im)
    // b1 = (a0.re - a1.re + a2.im - a3.im) + i (a0.im - a1.im - a2.re + a3.re)
    // b2 = (a0.re + a1.re - a2.re - a3.re) + i (a0.im + a1.im - a2.im - a3.im)
    // b3 = (a0.re - a1.re - a2.im + a3.im) + i (a0.im - a1.im + a2.re - a3.re)
    //
    // b0 = a0 + a1 +   a2 +   a3
    // b1 = a0 - a1 - i a2 + i a3
    // b2 = a0 + a2 -   a2 -   a3
    // b3 = a0 - a1 + i a2 - i a3

    let mut t1 = a[1].re;
    let mut t3 = a[0].re - t1;
    let mut t6 = a[2].re;
    t1 += a[0].re;
    let t8 = a[3].re - t6;
    t6 += a[3].re;
    a[0].re = t1 + t6;
    t1 -= t6;
    a[2].re = t1;

    let mut t2 = a[1].im;
    let mut t4 = a[0].im - t2;
    t2 += a[0].im;
    let mut t5 = a[3].im;
    a[1].im = t4 + t8;
    t4 -= t8;
    a[3].im = t4;

    let t7 = a[2].im - t5;
    t5 += a[2].im;
    a[1].re = t3 + t7;
    t3 -= t7;
    a[3].re = t3;
    a[0].im = t2 + t5;
    t2 -= t5;
    a[2].im = t2;

    normalise_all(a);
}

pub(crate) fn u8(a: &mut [Complex]) {
    u4(&mut a[..4]);

    let mut t1 = a[5].re;
    a[5].re = a[4].re - t1;
    t1 += a[4].re;

    let mut t3 = a[7].re;
    a[7].re = a[6].re - t3;
    t3 += a[6].re;

    let mut t8 = t3 - t1;
    t1 += t3;

    let mut t6 = a[2].im - t8;
    t8 += a[2].im;
    a[2].im = t8;

    let mut t5 = a[0].re - t1;
    a[4].re = t5;
    t1 += a[0].re;
    a[0].re = t1;

    let mut t2 = a[5].im;
    a[5].im = a[4].im - t2;
    t2 += a[4].im;

    let mut t4 = a[7].im;
    a[7].im = a[6].im - t4;
    t4 += a[6].im;

    a[6].im = t6;

    let mut t7 = t2 - t4;
    t2 += t4;

    t3 = a[2].re - t7;
    a[6].re = t3;
    t7 += a[2].re;
    a[2].re = t7;

    t6 = a[0].im - t2;
    a[4].im = t6;
    t2 += a[0].im;
    a[0].im = t2;

    t6 = SQRTHALF;

    t1 = a[5].re;
    t2 = a[5].im - t1;
    t2 *= t6;
    t1 += a[5].im;
    t1 *= t6;
    t4 = a[7].im;
    t3 = a[7].re - t4;
    t3 *= t6;
    t4 += a[7].re;
    t4 *= t6;

    t8 = t3 - t1;
    t1 += t3;
    t7 = t2 - t4;
    t2 += t4;

    t4 = a[3].im - t8;
    a[7].im = t4;
    t5 = a[1].re - t1;
    a[5].re = t5;
    t3 = a[3].re - t7;
    a[7].re = t3;
    t6 = a[1].im - t2;
    a[5].im = t6;

    t8 += a[3].im;
    a[3].im = t8;
    t1 += a[1].re;
    a[1].re = t1;
    t7 += a[3].re;
    a[3].re = t7;
    t2 += a[1].im;
    a[1].im = t2;

    normalise_all(a);
}

pub(crate) fn u16(a: &mut [Complex]) {
    u8(&mut a[..8]);
    u4(&mut a[8..12]);
    u4(&mut a[12..]);

    // untransformzero(a[0], a[4], a[8], a[12]);
    // untransformhalf(a[2], a[6], a[10], a[14]);
    // untransform(a[1], a[5], a[9], a[13], d16[0].re, d16[0].im);
    // untransform(a[3], a[7], a[11], a[15], d16[0].im, d16[0].re);

    // TODO: This is some fugly shit...
    let mut a0 = a[0];
    let mut a1 = a[4];
    let mut a2 = a[8];
    let mut a3 = a[12];

    //untransformzero(&mut a[0], &mut a[4], &mut a[8], &mut a[12]);
    untransformzero(&mut a0, &mut a1, &mut a2, &mut a3);
    a[0] = a0;
    a[4] = a1;
    a[8] = a2;
    a[12] = a3;

    a0 = a[1];
    a1 = a[5];
    a2 = a[9];
    a3 = a[13];
    //untransform(&mut a[1], &mut a[5], &mut a[9], &mut a[13], D16[0].re, D16[0].im,);
    untransform(&mut a0, &mut a1, &mut a2, &mut a3, D16[0].re, D16[0].im);
    a[1] = a0;
    a[5] = a1;
    a[9] = a2;
    a[13] = a3;

    a0 = a[2];
    a1 = a[6];
    a2 = a[10];
    a3 = a[14];
    //untransformhalf(&mut a[2], &mut a[6], &mut a[10], &mut a[14]);
    untransformhalf(&mut a0, &mut a1, &mut a2, &mut a3);
    a[2] = a0;
    a[6] = a1;
    a[10] = a2;
    a[14] = a3;

    a0 = a[3];
    a1 = a[7];
    a2 = a[11];
    a3 = a[15];
    //untransform(&mut a[3], &mut a[7], &mut a[11], &mut a[15], D16[0].im, D16[0].re,);
    untransform(&mut a0, &mut a1, &mut a2, &mut a3, D16[0].im, D16[0].re);
    a[3] = a0;
    a[7] = a1;
    a[11] = a2;
    a[15] = a3;
}

/* a[0...8n-1], w[0...2n-2] */
fn upass(a: &mut [Complex], w: &[Complex]) {
    debug_assert_eq!(a.len() % 8, 0);

    let n = a.len() / 8;

    debug_assert!(n >= 2);
    debug_assert_eq!(w.len(), 2 * n - 1);

    // Split a into four chunks of size 2*n.
    let (a, a1) = a.split_at_mut(2 * n);
    let (a1, a2) = a1.split_at_mut(2 * n);
    let (a2, a3) = a2.split_at_mut(2 * n);

    untransformzero(&mut a[0], &mut a1[0], &mut a2[0], &mut a3[0]);

    // TODO: Can I not use transformhalf here for some i?

    // TODO: The original version pulled the first iteration out of
    // the loop and unrolled the loop two iterations; check whether
    // that actually improves things here.
    for i in 1..2 * n {
        untransform(
            &mut a[i],
            &mut a1[i],
            &mut a2[i],
            &mut a3[i],
            w[i - 1].re,
            w[i - 1].im,
        );
    }
}

pub(crate) fn backward32(a: &mut [Complex]) {
    debug_assert_eq!(a.len(), 32);

    u16(&mut a[..16]);
    u8(&mut a[16..24]);
    u8(&mut a[24..32]);
    upass(a, &D32);
}

pub(crate) fn backward64(a: &mut [Complex]) {
    debug_assert_eq!(a.len(), 64);

    backward32(&mut a[..32]);
    u16(&mut a[32..48]);
    u16(&mut a[48..64]);
    upass(a, &D64);
}

pub(crate) fn backward128(a: &mut [Complex]) {
    debug_assert_eq!(a.len(), 128);

    backward64(&mut a[..64]);
    backward32(&mut a[64..96]);
    backward32(&mut a[96..128]);
    upass(a, &D128);
}

pub(crate) fn backward256(a: &mut [Complex]) {
    debug_assert_eq!(a.len(), 256);

    backward128(&mut a[..128]);
    backward64(&mut a[128..192]);
    backward64(&mut a[192..256]);
    upass(a, &D256);
}

pub(crate) fn backward512(a: &mut [Complex]) {
    debug_assert_eq!(a.len(), 512);

    backward256(&mut a[..256]);
    backward128(&mut a[256..384]);
    backward128(&mut a[384..512]);
    upass(a, &D512);
}

pub(crate) fn backward1024(a: &mut [Complex]) {
    debug_assert_eq!(a.len(), 1024);

    backward512(&mut a[..512]);
    backward256(&mut a[512..768]);
    backward256(&mut a[768..1024]);
    upass(a, &D1024);
}

pub(crate) fn backward2048(a: &mut [Complex]) {
    debug_assert_eq!(a.len(), 2048);

    backward1024(&mut a[..1024]);
    backward512(&mut a[1024..1536]);
    backward512(&mut a[1536..2048]);
    upass(a, &D2048);
}

pub(crate) fn backward4096(a: &mut [Complex]) {
    debug_assert_eq!(a.len(), 4096);

    backward2048(&mut a[..2048]);
    backward1024(&mut a[2048..3072]);
    backward1024(&mut a[3072..4096]);
    upass(a, &D4096);
}
