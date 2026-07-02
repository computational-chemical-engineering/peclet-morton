/// @file hd.hpp
/// @brief The `MORTON_HD` host/device function marker, shared by every core header.
///
/// Defining the marker in one tiny header (rather than inline in `morton.hpp`) lets headers that
/// `morton.hpp` includes *before* its class body — notably `wide_uint.hpp` — annotate their own
/// operators as device-callable too. Without this, `wide_uint`'s arithmetic could not run inside a
/// Kokkos kernel, which would break a device-resident Morton path for codes wider than 128 bits.
///
/// When the Morton Kokkos backend is enabled and Kokkos is already in the translation unit, defer
/// to Kokkos's own function marker: `KOKKOS_FUNCTION` resolves to the right host/device attributes
/// for whichever backend (CUDA / HIP / OpenMP / Serial) Kokkos was built with — so the same code is
/// callable from a Kokkos kernel on any backend. The raw `__CUDACC__`/`__HIPCC__` branch remains
/// the fallback for a direct nvcc/hipcc compile without Kokkos; an ordinary host build expands it
/// to nothing, leaving the CPU library byte-for-byte unchanged.
///
/// SPDX-License-Identifier: MIT

#ifndef MORTON_HD_HPP
#define MORTON_HD_HPP

#if !defined(MORTON_HD)
#if defined(MORTON_ENABLE_KOKKOS) && defined(KOKKOS_VERSION)
#define MORTON_HD KOKKOS_FUNCTION
#elif defined(__CUDACC__) || defined(__HIPCC__)
#define MORTON_HD __host__ __device__
#else
#define MORTON_HD
#endif
#endif  // !defined(MORTON_HD)

#endif  // MORTON_HD_HPP
