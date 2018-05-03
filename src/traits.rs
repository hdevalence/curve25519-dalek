// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// Copyright (c) 2016-2018 Isis Lovecruft, Henry de Valence
// See LICENSE for licensing information.
//
// Authors:
// - Isis Agora Lovecruft <isis@patternsinthevoid.net>
// - Henry de Valence <hdevalence@hdevalence.ca>

//! Module for common traits.

use core::borrow::Borrow;

use subtle;

use scalar::Scalar;

// ------------------------------------------------------------------------
// Public Traits
// ------------------------------------------------------------------------

/// Trait for getting the identity element of a point type.
pub trait Identity {
    /// Returns the identity element of the curve.
    /// Can be used as a constructor.
    fn identity() -> Self;
}

/// Trait for testing if a curve point is equivalent to the identity point.
pub trait IsIdentity {
    /// Return true if this element is the identity element of the curve.
    fn is_identity(&self) -> bool;
}

/// Implement generic identity equality testing for a point representations
/// which have constant-time equality testing and a defined identity
/// constructor.
impl<T> IsIdentity for T
where
    T: subtle::ConstantTimeEq + Identity,
{
    fn is_identity(&self) -> bool {
        self.ct_eq(&T::identity()).unwrap_u8() == 1u8
    }
}

/// A trait for constant-time multiscalar multiplication without precomputation.
pub trait MultiscalarMul {
    /// The type of point being multiplied, e.g., `RistrettoPoint`.
    type Point;

    /// Given an iterator of (possibly secret) scalars and an iterator of
    /// public points, compute
    /// $$
    /// Q = c\_1 P\_1 + \cdots + c\_n P\_n.
    /// $$
    ///
    /// It is an error to call this function with two iterators of different lengths.
    ///
    /// # Examples
    ///
    /// The trait bound aims for maximum flexibility: the inputs must be
    /// convertable to iterators (`I: IntoIter`), and the iterator's items
    /// must be `Borrow<Scalar>` (or `Borrow<Point>`), to allow
    /// iterators returning either `Scalar`s or `&Scalar`s.
    ///
    /// ```
    /// use curve25519_dalek::constants;
    /// use curve25519_dalek::traits::MultiscalarMul;
    /// use curve25519_dalek::ristretto::RistrettoPoint;
    /// use curve25519_dalek::scalar::Scalar;
    ///
    /// // Some scalars
    /// let a = Scalar::from_u64(87329482);
    /// let b = Scalar::from_u64(37264829);
    /// let c = Scalar::from_u64(98098098);
    ///
    /// // Some points
    /// let P = constants::RISTRETTO_BASEPOINT_POINT;
    /// let Q = P + P;
    /// let R = P + Q;
    ///
    /// // A1 = a*P + b*Q + c*R
    /// let abc = [a,b,c];
    /// let A1 = RistrettoPoint::multiscalar_mul(&abc, &[P,Q,R]);
    /// // Note: (&abc).into_iter(): Iterator<Item=&Scalar>
    ///
    /// // A2 = (-a)*P + (-b)*Q + (-c)*R
    /// let minus_abc = abc.iter().map(|x| -x);
    /// let A2 = RistrettoPoint::multiscalar_mul(minus_abc, &[P,Q,R]);
    /// // Note: minus_abc.into_iter(): Iterator<Item=Scalar>
    ///
    /// assert_eq!(A1.compress(), (-A2).compress());
    /// ```
    fn multiscalar_mul<I, J>(scalars: I, points: J) -> Self::Point
    where
        I: IntoIterator,
        I::Item: Borrow<Scalar>,
        J: IntoIterator,
        J::Item: Borrow<Self::Point>;
}

/// A trait for variable-time multiscalar multiplication without precomputation.
pub trait VartimeMultiscalarMul {
    /// The type of point being multiplied, e.g., `RistrettoPoint`.
    type Point;

    /// Given an iterator of (possibly secret) scalars and an iterator of
    /// public points, compute
    /// $$
    /// Q = c\_1 P\_1 + \cdots + c\_n P\_n.
    /// $$
    ///
    /// It is an error to call this function with two iterators of different lengths.
    ///
    /// # Examples
    ///
    /// The trait bound aims for maximum flexibility: the inputs must be
    /// convertable to iterators (`I: IntoIter`), and the iterator's items
    /// must be `Borrow<Scalar>` (or `Borrow<Point>`), to allow
    /// iterators returning either `Scalar`s or `&Scalar`s.
    ///
    /// ```
    /// use curve25519_dalek::constants;
    /// use curve25519_dalek::traits::MultiscalarMul;
    /// use curve25519_dalek::ristretto::RistrettoPoint;
    /// use curve25519_dalek::scalar::Scalar;
    ///
    /// // Some scalars
    /// let a = Scalar::from_u64(87329482);
    /// let b = Scalar::from_u64(37264829);
    /// let c = Scalar::from_u64(98098098);
    ///
    /// // Some points
    /// let P = constants::RISTRETTO_BASEPOINT_POINT;
    /// let Q = P + P;
    /// let R = P + Q;
    ///
    /// // A1 = a*P + b*Q + c*R
    /// let abc = [a,b,c];
    /// let A1 = RistrettoPoint::multiscalar_mul(&abc, &[P,Q,R]);
    /// // Note: (&abc).into_iter(): Iterator<Item=&Scalar>
    ///
    /// // A2 = (-a)*P + (-b)*Q + (-c)*R
    /// let minus_abc = abc.iter().map(|x| -x);
    /// let A2 = RistrettoPoint::multiscalar_mul(minus_abc, &[P,Q,R]);
    /// // Note: minus_abc.into_iter(): Iterator<Item=Scalar>
    ///
    /// assert_eq!(A1.compress(), (-A2).compress());
    /// ```
    fn vartime_multiscalar_mul<I, J>(scalars: I, points: J) -> Self::Point
    where
        I: IntoIterator,
        I::Item: Borrow<Scalar>,
        J: IntoIterator,
        J::Item: Borrow<Self::Point>;
}

/// A trait for constant-time multiscalar multiplication with precomputation.
///
/// A general multiscalar multiplication with precomputation can be written as
/// $$
/// Q = a_1 A_1 + \cdots + a_n A_n + b_1 B_1 + \cdots + b_m B_m,
/// $$
/// where the \\(B_i\\) are *static* points, for which precomputation
/// is possible, and the \\(A_j\\) are *dynamic* points, for which
/// precomputation is not possible.
pub trait PrecomputedMultiscalarMul: Sized {
    /// The type of point to be multiplied, e.g., `RistrettoPoint`.
    type Point;

    /// Given the static points \\( B_i \\), perform precomputation
    /// and return the precomputation data.
    fn new<I>(static_points: I) -> Self
    where
        I: IntoIterator,
        I::Item: Borrow<Self::Point>;

    /// Given `static_scalars`, an iterator of (possibly secret)
    /// scalars \\(b_i\\), `dynamic_scalars`, an iterator of (possibly
    /// secret) scalars \\(a_i\\), and `dynamic_points`, an iterator
    /// of points \\(A_i\\), compute
    /// $$
    /// Q = a_1 A_1 + \cdots + a_n A_n + b_1 B_1 + \cdots + b_m B_m,
    /// $$
    /// where the \\(B_j\\) are the points that were supplied to `new`.
    ///
    /// It is an error to call this function with iterators of
    /// inconsistent lengths.
    ///
    /// The trait bound aims for maximum flexibility: the inputs must be
    /// convertable to iterators (`I: IntoIter`), and the iterator's items
    /// must be `Borrow<Scalar>` (or `Borrow<Point>`), to allow
    /// iterators returning either `Scalar`s or `&Scalar`s.
    fn mixed_multiscalar_mul<I, J, K>(
        &self,
        static_scalars: I,
        dynamic_scalars: J,
        dynamic_points: K,
    ) -> Self::Point
    where
        I: IntoIterator,
        I::Item: Borrow<Scalar>,
        J: IntoIterator,
        J::Item: Borrow<Scalar>,
        K: IntoIterator,
        K::Item: Borrow<Self::Point>;

    /// Given `static_scalars`, an iterator of (possibly secret)
    /// scalars \\(b_i\\), compute
    /// $$
    /// Q = b_1 B_1 + \cdots + b_m B_m,
    /// $$
    /// where the \\(B_j\\) are the points that were supplied to `new`.
    ///
    /// It is an error to call this function with iterators of
    /// inconsistent lengths.
    ///
    /// The trait bound aims for maximum flexibility: the input must
    /// be convertable to iterators (`I: IntoIter`), and the
    /// iterator's items must be `Borrow<Scalar>`, to allow iterators
    /// returning either `Scalar`s or `&Scalar`s.
    fn multiscalar_mul<I>(&self, static_scalars: I) -> Self::Point
    where
        I: IntoIterator,
        I::Item: Borrow<Scalar>,
    {
        use std::iter;

        Self::mixed_multiscalar_mul(
            self,
            static_scalars,
            iter::empty::<Scalar>(),
            iter::empty::<Self::Point>(),
        )
    }
}

/// A trait for variable-time multiscalar multiplication with precomputation.
///
/// A general multiscalar multiplication with precomputation can be written as
/// $$
/// Q = a_1 A_1 + \cdots + a_n A_n + b_1 B_1 + \cdots + b_m B_m,
/// $$
/// where the \\(B_i\\) are *static* points, for which precomputation
/// is possible, and the \\(A_j\\) are *dynamic* points, for which
/// precomputation is not possible.
pub trait VartimePrecomputedMultiscalarMul: Sized {
    /// The type of point to be multiplied, e.g., `RistrettoPoint`.
    type Point;

    /// Given the static points \\( B_i \\), perform precomputation
    /// and return the precomputation data.
    fn new<I>(static_points: I) -> Self
    where
        I: IntoIterator,
        I::Item: Borrow<Self::Point>;

    /// Given `static_scalars`, an iterator of public scalars
    /// \\(b_i\\), `dynamic_scalars`, an iterator of public scalars
    /// \\(a_i\\), and `dynamic_points`, an iterator of points
    /// \\(A_i\\), compute
    /// $$
    /// Q = a_1 A_1 + \cdots + a_n A_n + b_1 B_1 + \cdots + b_m B_m,
    /// $$
    /// where the \\(B_j\\) are the points that were supplied to `new`.
    ///
    /// It is an error to call this function with iterators of
    /// inconsistent lengths.
    ///
    /// The trait bound aims for maximum flexibility: the inputs must be
    /// convertable to iterators (`I: IntoIter`), and the iterator's items
    /// must be `Borrow<Scalar>` (or `Borrow<Point>`), to allow
    /// iterators returning either `Scalar`s or `&Scalar`s.
    fn vartime_mixed_multiscalar_mul<I, J, K>(
        &self,
        static_scalars: I,
        dynamic_scalars: J,
        dynamic_points: K,
    ) -> Self::Point
    where
        I: IntoIterator,
        I::Item: Borrow<Scalar>,
        J: IntoIterator,
        J::Item: Borrow<Scalar>,
        K: IntoIterator,
        K::Item: Borrow<Self::Point>;

    /// Given `static_scalars`, an iterator of public scalars
    /// \\(b_i\\), compute
    /// $$
    /// Q = b_1 B_1 + \cdots + b_m B_m,
    /// $$
    /// where the \\(B_j\\) are the points that were supplied to `new`.
    ///
    /// It is an error to call this function with iterators of
    /// inconsistent lengths.
    ///
    /// The trait bound aims for maximum flexibility: the input must
    /// be convertable to iterators (`I: IntoIter`), and the
    /// iterator's items must be `Borrow<Scalar>`, to allow iterators
    /// returning either `Scalar`s or `&Scalar`s.
    fn vartime_multiscalar_mul<I>(&self, static_scalars: I) -> Self::Point
    where
        I: IntoIterator,
        I::Item: Borrow<Scalar>,
    {
        use std::iter;

        Self::vartime_mixed_multiscalar_mul(
            self,
            static_scalars,
            iter::empty::<Scalar>(),
            iter::empty::<Self::Point>(),
        )
    }
}

// ------------------------------------------------------------------------
// Private Traits
// ------------------------------------------------------------------------

/// Trait for checking whether a point is on the curve.
///
/// This trait is only for debugging/testing, since it should be
/// impossible for a `curve25519-dalek` user to construct an invalid
/// point.
pub(crate) trait ValidityCheck {
    /// Checks whether the point is on the curve. Not CT.
    fn is_valid(&self) -> bool;
}
