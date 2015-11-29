## radixsort-rs

[![Build Status](https://travis-ci.org/bitshifter/radixsort-rs.svg?branch=master)](https://travis-ci.org/bitshifter/radixsort-rs)

This is an implementation of radix sort in Rust.

**NOTE** the test and benchmark code is in the process of being moved from a
common format used by my C and C++ versions of radixsort to a more standard
Rust setup.

The radix sort implentation also include a 11 bit radix sort which uses less
passes than the standard 8 bit radix sort, and a floating point key radix sort,
as described by Michael Herf:

 * http://stereopsis.com/radix.html

The radix sort interface is designed to perform no memory allocations, and does
not write the final sorted values to an output buffer, but rather returns which
buffer the final result resides in. The rationaly here is the callig code can
decide if it's appropriate to allocate temporary space, or to copy the results
to a specified buffer. If so desired a simple wrapper could be written to
simplify the interface at the expense of some performance and memory
allocation.

## License

This software is licensed under the zlib license, see the LICENSE file for
details.
