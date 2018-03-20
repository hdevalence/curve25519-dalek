#![allow(non_snake_case)]

use core::ops::{Add, Sub};

use scalar::Scalar;
use edwards::EdwardsPoint;
use curve_models::CompletedPoint;

use subtle::ConditionallyAssignable;
use subtle::ConditionallyNegatable;

use super::window::LookupTable;
use traits::{Doubleable, Identity};

/// Perform constant-time, variable-base scalar multiplication.
///
/// Type variables:
/// - `AddLHS`: the main point type
/// - `AddRHS`: a cached point type, used by the lookup table
/// - `Out`: the output of addition.
pub(crate) fn mul<AddLHS, AddRHS, Out>(point: &AddLHS, scalar: &Scalar) -> AddLHS
where
    for<'a, 'b> &'a AddLHS: Add<&'b AddRHS, Output = Out>,
    for<'a> LookupTable<AddRHS>: From<&'a AddLHS>,
    AddLHS: Identity + Doubleable + From<Out>,
    AddRHS: Identity + ConditionallyAssignable + ConditionallyNegatable,
{
    // Construct a lookup table of [P,2P,3P,4P,5P,6P,7P,8P]
    let lookup_table = LookupTable::<AddRHS>::from(point);
    // Setting s = scalar, compute
    //
    //    s = s_0 + s_1*16^1 + ... + s_63*16^63,
    //
    // with `-8 ≤ s_i < 8` for `0 ≤ i < 63` and `-8 ≤ s_63 ≤ 8`.
    let scalar_digits = scalar.to_radix_16();
    // Compute s*P as
    //
    //    s*P = P*(s_0 +   s_1*16^1 +   s_2*16^2 + ... +   s_63*16^63)
    //    s*P =  P*s_0 + P*s_1*16^1 + P*s_2*16^2 + ... + P*s_63*16^63
    //    s*P = P*s_0 + 16*(P*s_1 + 16*(P*s_2 + 16*( ... + P*s_63)...))
    //
    // We sum right-to-left.
    let mut Q = AddLHS::identity();
    for i in (0..64).rev() {
        // Q <-- 16*Q
        Q = Q.mul_by_pow_2(4);
        // Q <-- Q + P * s_i
        Q = (&Q + &lookup_table.select(scalar_digits[i])).into();
    }
    Q
}
