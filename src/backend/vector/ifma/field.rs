// -*- mode: rust; coding: utf-8; -*-
//
// This file is part of curve25519-dalek.
// Copyright (c) 2018 Henry de Valence
// See LICENSE for licensing information.
//
// Authors:
// - Henry de Valence <hdevalence@hdevalence.ca>

use core::ops::{Add, Mul, Neg};
use packed_simd::{i32x8, u32x8, u64x4, IntoBits};

use backend::u64::field::FieldElement51;

