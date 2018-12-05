An AVX512-IFMA implementation of the vectorized point operation
strategy.

# IFMA instructions

## IFMA for big-integer multiplications


## An alternative strategy

The strategy described above is aimed at big-integer multiplications,
such as 1024, 2048, or 4096 bits, which would be used for applications
like RSA.  However, elliptic curve cryptography uses much smaller field
sizes, such as 256 or 384 bits, so a different strategy is needed.

Instead, because we have parallelism at the level of the formulas for
curve operations, we can vectorize *across* field elements, using part
of each vector for a different field element.  The parallel Edwards
formulas provide 4-way parallelism, so they can be implemented using
256-bit vectors using a single 64-bit lane for each element, or using
512-bit vectors using two 64-bit lanes.

This means that inside the field computation, it's only necessary to
achieve 2-way parallelism, and the data for each field element can be
contained in 128-bit lanes, so that cross-lane operations can use the
faster `vpshufd` (1c latency) instead of a general shuffle instruction
(3c latency). However, the only available CPU supporting IFMA (the
i3-8121U) executes 512-bit IFMA instructions at half rate compared to
256-bit instructions, so for now there's no throughput advantage to
using 512-bit IFMA instructions, and this implementation uses 256-bit
vectors.

Using this approach, instead of scanning through the terms of the
source operands as above, we arrange the computation in
product-scanning form, as described below.

# Multiplication

# Choice of radix

# Reductions


# Squaring

[2016_gueron_krasnov]: https://ieeexplore.ieee.org/document/7563269
[2018_drucker_gueron]: https://eprint.iacr.org/2018/335
