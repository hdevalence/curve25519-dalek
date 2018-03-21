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

use core::ops::Add;
use core::borrow::Borrow;

use clear_on_drop::ClearOnDrop;

use subtle::ConditionallyAssignable;
use subtle::ConditionallyNegatable;

use scalar::Scalar;
use edwards::EdwardsPoint;
use curve_models::CompletedPoint;

use super::window::LookupTable;
use traits::{Doubleable, Identity};

/// Perform constant-time, variable-base scalar multiplication.
///
/// Type variables:
/// - `AddLHS`: the main point type
/// - `AddRHS`: a cached point type, used by the lookup table
/// - `Out`: the output of addition.
pub(crate) fn multiscalar_mul<I, J, AddLHS, AddRHS, Out>(scalars: I, points: J) -> AddLHS
where
    I: IntoIterator,
    I::Item: Borrow<Scalar>,
    J: IntoIterator,
    J::Item: Borrow<AddLHS>,
    for<'a, 'b> &'a AddLHS: Add<&'b AddRHS, Output = Out>,
    for<'a> LookupTable<AddRHS>: From<&'a AddLHS>,
    AddLHS: Identity + Doubleable + From<Out>,
    AddRHS: Identity + ConditionallyAssignable + ConditionallyNegatable,
{
    // Construct a lookup table of [P,2P,3P,4P,5P,6P,7P,8P]
    // for each input point P
    let lookup_tables: Vec<_> = points
        .into_iter()
        .map(|point| LookupTable::<AddRHS>::from(point.borrow()))
        .collect();

    // Setting s_i = i-th scalar, compute
    //
    //    s_i = s_{i,0} + s_{i,1}*16^1 + ... + s_{i,63}*16^63,
    //
    // with `-8 ≤ s_{i,j} < 8` for `0 ≤ j < 63` and `-8 ≤ s_{i,63} ≤ 8`.
    //
    // This puts the scalar digits into a heap-allocated Vec.
    // To ensure that these are erased, pass ownership of the Vec into a
    // ClearOnDrop wrapper.
    let scalar_digits_vec: Vec<_> = scalars
        .into_iter()
        .map(|s| s.borrow().to_radix_16())
        .collect();
    let scalar_digits = ClearOnDrop::new(scalar_digits_vec);

    // Compute s_1*P_1 + ... + s_n*P_n: since
    //
    //    s_i*P_i = P_i*(s_{i,0} +     s_{i,1}*16^1 + ... +     s_{i,63}*16^63)
    //    s_i*P_i =  P_i*s_{i,0} + P_i*s_{i,1}*16^1 + ... + P_i*s_{i,63}*16^63
    //    s_i*P_i =  P_i*s_{i,0} + 16*(P_i*s_{i,1} + 16*( ... + 16*P_i*s_{i,63})...)
    //
    // we have the two-dimensional sum
    //
    //    s_1*P_1 =   P_1*s_{1,0} + 16*(P_1*s_{1,1} + 16*( ... + 16*P_1*s_{1,63})...)
    //  + s_2*P_2 = + P_2*s_{2,0} + 16*(P_2*s_{2,1} + 16*( ... + 16*P_2*s_{2,63})...)
    //      ...
    //  + s_n*P_n = + P_n*s_{n,0} + 16*(P_n*s_{n,1} + 16*( ... + 16*P_n*s_{n,63})...)
    //
    // We sum column-wise top-to-bottom, then right-to-left,
    // multiplying by 16 only once per column.
    //
    // This provides the speedup over doing n independent scalar
    // mults: we perform 63 multiplications by 16 instead of 63*n
    // multiplications, saving 252*(n-1) doublings.
    let mut Q = AddLHS::identity();
    for j in (0..64).rev() {
        Q = Q.mul_by_pow_2(4);
        let it = scalar_digits.iter().zip(lookup_tables.iter());
        for (s_i, lookup_table_i) in it {
            // R_i = s_{i,j} * P_i
            let R_i = lookup_table_i.select(s_i[j]);
            // Q = Q + R_i
            Q = (&Q + &R_i).into();
        }
    }
    Q
}
