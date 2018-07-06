// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// Copyright (c) 2016-2018 Isis Lovecruft, Henry de Valence
// See LICENSE for licensing information.
//
// Authors:
// - Isis Agora Lovecruft <isis@patternsinthevoid.net>
// - Henry de Valence <hdevalence@hdevalence.ca>

#![allow(non_snake_case)]

use core::borrow::Borrow;

use clear_on_drop::ClearOnDrop;

use backend::avx2::edwards::{CachedPoint, ExtendedPoint};
use edwards::EdwardsPoint;
use scalar::Scalar;
use scalar_mul::window::{LookupTable, NafLookupTable5, NafLookupTable8};
use traits::Identity;
use traits::MultiscalarMul;
use traits::PrecomputedMultiscalarMul;
use traits::VartimeMultiscalarMul;
use traits::VartimePrecomputedMultiscalarMul;

/// Multiscalar multiplication using interleaved window / Straus'
/// method.  See the `Straus` struct in the serial backend for more
/// details.
///
/// This exists as a seperate implementation from that one because the
/// AVX2 code uses different curve models (it does not pass between
/// multiple models during scalar mul), and it has to convert the
/// point representation on the fly.
pub struct Straus {}

#[cfg(any(feature = "alloc", feature = "std"))]
impl MultiscalarMul for Straus {
    type Point = EdwardsPoint;

    fn multiscalar_mul<I, J>(scalars: I, points: J) -> EdwardsPoint
    where
        I: IntoIterator,
        I::Item: Borrow<Scalar>,
        J: IntoIterator,
        J::Item: Borrow<EdwardsPoint>,
    {
        // Construct a lookup table of [P,2P,3P,4P,5P,6P,7P,8P]
        // for each input point P
        let lookup_tables: Vec<_> = points
            .into_iter()
            .map(|point| LookupTable::<CachedPoint>::from(point.borrow()))
            .collect();

        let scalar_digits_vec: Vec<_> = scalars
            .into_iter()
            .map(|s| s.borrow().to_radix_16())
            .collect();
        // Pass ownership to a ClearOnDrop wrapper
        let scalar_digits = ClearOnDrop::new(scalar_digits_vec);

        let mut Q = ExtendedPoint::identity();
        for j in (0..64).rev() {
            Q = Q.mul_by_pow_2(4);
            let it = scalar_digits.iter().zip(lookup_tables.iter());
            for (s_i, lookup_table_i) in it {
                // Q = Q + s_{i,j} * P_i
                Q = &Q + &lookup_table_i.select(s_i[j]);
            }
        }
        Q.into()
    }
}

#[cfg(any(feature = "alloc", feature = "std"))]
impl VartimeMultiscalarMul for Straus {
    type Point = EdwardsPoint;

    fn vartime_multiscalar_mul<I, J>(scalars: I, points: J) -> EdwardsPoint
    where
        I: IntoIterator,
        I::Item: Borrow<Scalar>,
        J: IntoIterator,
        J::Item: Borrow<EdwardsPoint>,
    {
        let nafs: Vec<_> = scalars
            .into_iter()
            .map(|c| c.borrow().non_adjacent_form(5))
            .collect();
        let lookup_tables: Vec<_> = points
            .into_iter()
            .map(|point| NafLookupTable5::<CachedPoint>::from(point.borrow()))
            .collect();

        let mut Q = ExtendedPoint::identity();

        for i in (0..255).rev() {
            Q = Q.double();

            for (naf, lookup_table) in nafs.iter().zip(lookup_tables.iter()) {
                if naf[i] > 0 {
                    Q = &Q + &lookup_table.select(naf[i] as usize);
                } else if naf[i] < 0 {
                    Q = &Q - &lookup_table.select(-naf[i] as usize);
                }
            }
        }
        Q.into()
    }
}

/// Precomputation for constant-time multiscalar multiplication using Straus' algorithm.
#[cfg(any(feature = "alloc", feature = "std"))]
pub struct PrecomputedStraus {
    lookup_tables: Vec<LookupTable<CachedPoint>>,
}

#[cfg(any(feature = "alloc", feature = "std"))]
impl PrecomputedMultiscalarMul for PrecomputedStraus {
    type Point = EdwardsPoint;

    fn new<I>(static_points: I) -> Self
    where
        I: IntoIterator,
        I::Item: Borrow<Self::Point>,
    {
        PrecomputedStraus {
            lookup_tables: static_points
                .into_iter()
                .map(|point| LookupTable::<CachedPoint>::from(point.borrow()))
                .collect(),
        }
    }

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
        K::Item: Borrow<Self::Point>,
    {
        // Compute the scalar digits for the static and dynamic scalars.
        // To ensure that these are erased, pass ownership of the Vec into a
        // ClearOnDrop wrapper.
        use clear_on_drop::ClearOnDrop;

        let static_scalar_digits_vec: Vec<_> = static_scalars
            .into_iter()
            .map(|s| s.borrow().to_radix_16())
            .collect();
        let static_scalar_digits = ClearOnDrop::new(static_scalar_digits_vec);

        let dynamic_scalar_digits_vec: Vec<_> = dynamic_scalars
            .into_iter()
            .map(|s| s.borrow().to_radix_16())
            .collect();
        let dynamic_scalar_digits = ClearOnDrop::new(dynamic_scalar_digits_vec);

        let dynamic_lookup_tables: Vec<_> = dynamic_points
            .into_iter()
            .map(|point| LookupTable::<CachedPoint>::from(point.borrow()))
            .collect();

        let mut Q = ExtendedPoint::identity();
        for j in (0..64).rev() {
            Q = Q.mul_by_pow_2(4);

            // Add the static points
            let it = static_scalar_digits.iter().zip(self.lookup_tables.iter());
            for (s_i, lookup_table_i) in it {
                // Q = Q + s_{i,j} * P_i
                Q = &Q + &lookup_table_i.select(s_i[j]);
            }

            // Add the dynamic points
            let it = dynamic_scalar_digits
                .iter()
                .zip(dynamic_lookup_tables.iter());
            for (s_i, lookup_table_i) in it {
                // Q = Q + s_{i,j} * P_i
                Q = &Q + &lookup_table_i.select(s_i[j]);
            }
        }

        Q.into()
    }
}

/// Precomputation for variable-time multiscalar multiplication using Straus' algorithm.
#[cfg(any(feature = "alloc", feature = "std"))]
pub struct VartimePrecomputedStraus {
    static_tables: Vec<NafLookupTable8<CachedPoint>>,
}

#[cfg(any(feature = "alloc", feature = "std"))]
impl VartimePrecomputedMultiscalarMul for VartimePrecomputedStraus {
    type Point = EdwardsPoint;

    fn new<I>(static_points: I) -> Self
    where
        I: IntoIterator,
        I::Item: Borrow<Self::Point>,
    {
        VartimePrecomputedStraus {
            static_tables: static_points
                .into_iter()
                .map(|point| NafLookupTable8::<CachedPoint>::from(point.borrow()))
                .collect(),
        }
    }

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
        K::Item: Borrow<Self::Point>,
    {
        // Scalars for precomputed points use a NAF of width 8
        let static_nafs: Vec<_> = static_scalars
            .into_iter()
            .map(|s| s.borrow().non_adjacent_form(8))
            .collect();

        // Scalars for dynamic points use a NAF of width 5
        let dynamic_nafs: Vec<_> = dynamic_scalars
            .into_iter()
            .map(|s| s.borrow().non_adjacent_form(5))
            .collect();

        // Build lookup tables for dynamic points
        let dynamic_tables: Vec<_> = dynamic_points
            .into_iter()
            .map(|point| NafLookupTable5::<CachedPoint>::from(point.borrow()))
            .collect();

        let mut Q = ExtendedPoint::identity();

        for i in (0..255).rev() {
            Q = Q.double();

            for (naf, table) in static_nafs.iter().zip(self.static_tables.iter()) {
                if naf[i] > 0 {
                    Q = &Q + &table.select(naf[i] as usize);
                } else if naf[i] < 0 {
                    Q = &Q - &table.select(-naf[i] as usize);
                }
            }

            for (naf, table) in dynamic_nafs.iter().zip(dynamic_tables.iter()) {
                if naf[i] > 0 {
                    Q = &Q + &table.select(naf[i] as usize);
                } else if naf[i] < 0 {
                    Q = &Q - &table.select(-naf[i] as usize);
                }
            }
        }

        Q.into()
    }
}

#[cfg(all(test, feature = "stage2_build"))]
mod test {
    use super::*;

    use rand::OsRng;

    use constants;
    use edwards::EdwardsPoint;

    #[test]
    fn multiscalar_mul_consistency() {
        let n = 32;

        let mut rng = OsRng::new().unwrap();

        let static_scalars: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();
        let dynamic_scalars: Vec<Scalar> = (0..n).map(|_| Scalar::random(&mut rng)).collect();

        let static_points: Vec<EdwardsPoint> = static_scalars
            .iter()
            .map(|s| s * &constants::ED25519_BASEPOINT_TABLE)
            .collect();
        let dynamic_points: Vec<EdwardsPoint> = dynamic_scalars
            .iter()
            .map(|s| s * &constants::ED25519_BASEPOINT_TABLE)
            .collect();

        let ct_precomp = PrecomputedStraus::new(&static_points);
        let vt_precomp = VartimePrecomputedStraus::new(&static_points);

        let res1 = EdwardsPoint::multiscalar_mul(
            static_scalars.iter().chain(dynamic_scalars.iter()),
            static_points.iter().chain(dynamic_points.iter()),
        );

        let res2 =
            ct_precomp.mixed_multiscalar_mul(&static_scalars, &dynamic_scalars, &dynamic_points);

        let res3 = vt_precomp.vartime_mixed_multiscalar_mul(
            &static_scalars,
            &dynamic_scalars,
            &dynamic_points,
        );

        assert_eq!(res1, res2);
        assert_eq!(res1, res3);
    }
}
