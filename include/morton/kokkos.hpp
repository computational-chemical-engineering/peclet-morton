/// @file kokkos.hpp
/// @brief Portable Kokkos backend for the morton-arithmetic library.
///
/// The per-element work is the ordinary `Morton<Dim,Bits>` code path (marked
/// `MORTON_HD`, which resolves to `KOKKOS_FUNCTION` once Kokkos is in the
/// translation unit), so the device kernels share the exact, already-tested
/// logic of the CPU library -- there is no second implementation to keep in
/// sync. Unlike the retired raw-CUDA backend, this runs on **any** Kokkos
/// backend (CUDA / HIP / OpenMP / Serial); the backend and architecture come
/// from the Kokkos build the suite is pointed at, not from this header.
///
/// Two layers, mirroring the rest of the suite's Kokkos idiom:
///   * View-based ops take caller-owned `Kokkos::View`s -- these compose into
///     larger pipelines and are exactly what the transport-core halo machinery
///     moves between ranks. The execution space is deduced from the View.
///   * `*_host` convenience wrappers take raw host pointers + `n`; they stage
///     into device Views, run, and copy back. Handy for one-shot calls/tests.
///
/// Encode/decode on the device necessarily use the software (PDEP-free) bit path
/// (`deposit`/`extract` guard the x86 intrinsics behind `!__CUDA_ARCH__`).
///
/// SPDX-License-Identifier: MIT

#ifndef MORTON_KOKKOS_HPP
#define MORTON_KOKKOS_HPP

#if !defined(MORTON_ENABLE_KOKKOS)
#error "morton/kokkos.hpp requires MORTON_ENABLE_KOKKOS (configure with -DMORTON_ENABLE_KOKKOS=ON)"
#endif

#include <cstddef>
#include <Kokkos_Core.hpp>  // must precede morton.hpp so KOKKOS_VERSION drives MORTON_HD

#include "morton/morton.hpp"

namespace morton {
namespace kokkos {

namespace detail {
// Unmanaged host view over a raw caller pointer, for staging to/from device.
template <class T>
using host_ptr_view = Kokkos::View<T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
}  // namespace detail

// ---------------------------------------------------------------------------
// View-based API (caller owns the Views; execution space deduced from `out`)
// ---------------------------------------------------------------------------

/// Encode coordinate Views (2D) into a code View. `x`,`y`,`out` extents must match.
template <unsigned Bits, class ViewX, class ViewY, class ViewOut>
void encode2(const ViewX& x, const ViewY& y, const ViewOut& out) {
  using M = Morton<2, Bits>;
  using exec = typename ViewOut::execution_space;
  const std::size_t n = out.extent(0);
  Kokkos::parallel_for(
      "morton::kokkos::encode2", Kokkos::RangePolicy<exec>(0, n),
      KOKKOS_LAMBDA(const std::size_t i) { out(i) = M::encode(x(i), y(i)).code(); });
}

/// Encode coordinate Views (3D) into a code View.
template <unsigned Bits, class ViewX, class ViewY, class ViewZ, class ViewOut>
void encode3(const ViewX& x, const ViewY& y, const ViewZ& z, const ViewOut& out) {
  using M = Morton<3, Bits>;
  using exec = typename ViewOut::execution_space;
  const std::size_t n = out.extent(0);
  Kokkos::parallel_for(
      "morton::kokkos::encode3", Kokkos::RangePolicy<exec>(0, n),
      KOKKOS_LAMBDA(const std::size_t i) { out(i) = M::encode(x(i), y(i), z(i)).code(); });
}

/// Decode a code View back to coordinate Views (2D).
template <unsigned Bits, class ViewCode, class ViewX, class ViewY>
void decode2(const ViewCode& code, const ViewX& x, const ViewY& y) {
  using M = Morton<2, Bits>;
  using exec = typename ViewCode::execution_space;
  const std::size_t n = code.extent(0);
  Kokkos::parallel_for(
      "morton::kokkos::decode2", Kokkos::RangePolicy<exec>(0, n),
      KOKKOS_LAMBDA(const std::size_t i) {
        M m = M::from_code(code(i));
        x(i) = m.get(0);
        y(i) = m.get(1);
      });
}

/// Decode a code View back to coordinate Views (3D).
template <unsigned Bits, class ViewCode, class ViewX, class ViewY, class ViewZ>
void decode3(const ViewCode& code, const ViewX& x, const ViewY& y, const ViewZ& z) {
  using M = Morton<3, Bits>;
  using exec = typename ViewCode::execution_space;
  const std::size_t n = code.extent(0);
  Kokkos::parallel_for(
      "morton::kokkos::decode3", Kokkos::RangePolicy<exec>(0, n),
      KOKKOS_LAMBDA(const std::size_t i) {
        M m = M::from_code(code(i));
        x(i) = m.get(0);
        y(i) = m.get(1);
        z(i) = m.get(2);
      });
}

/// Add a signed `delta` to `axis` of every code (wraps mod 2^Bits). `in`/`out`
/// may be the same View.
template <unsigned Dim, unsigned Bits, class ViewIn, class ViewOut>
void add(const ViewIn& in, const ViewOut& out, unsigned axis, long long delta) {
  using M = Morton<Dim, Bits>;
  using C = typename M::coord_type;
  using exec = typename ViewOut::execution_space;
  const std::size_t n = out.extent(0);
  const bool neg = delta < 0;
  const C mag =
      neg ? C(static_cast<unsigned long long>(-delta)) : C(static_cast<unsigned long long>(delta));
  Kokkos::parallel_for(
      "morton::kokkos::add", Kokkos::RangePolicy<exec>(0, n), KOKKOS_LAMBDA(const std::size_t i) {
        M m = M::from_code(in(i));
        if (neg)
          m.sub(axis, mag);
        else
          m.add(axis, mag);
        out(i) = m.code();
      });
}

/// Subtract `k` from `axis` of every code (wraps). Convenience over `add`.
template <unsigned Dim, unsigned Bits, class ViewIn, class ViewOut>
void sub(const ViewIn& in, const ViewOut& out, unsigned axis,
         typename Morton<Dim, Bits>::coord_type k) {
  add<Dim, Bits>(in, out, axis, -static_cast<long long>(k));
}

/// Step every code one cell along ±`axis` (dir>=0 -> +1, else -1).
template <unsigned Dim, unsigned Bits, class ViewIn, class ViewOut>
void step(const ViewIn& in, const ViewOut& out, unsigned axis, int dir) {
  add<Dim, Bits>(in, out, axis, dir >= 0 ? 1 : -1);
}

// ---------------------------------------------------------------------------
// Host-pointer convenience API (stage to device Views, run, copy back)
// ---------------------------------------------------------------------------

template <unsigned Bits>
void encode2_host(const typename Morton<2, Bits>::coord_type* x,
                  const typename Morton<2, Bits>::coord_type* y,
                  typename Morton<2, Bits>::code_type* out, std::size_t n) {
  using C = typename Morton<2, Bits>::coord_type;
  using K = typename Morton<2, Bits>::code_type;
  Kokkos::View<C*> dx("morton.x", n), dy("morton.y", n);
  Kokkos::View<K*> dout("morton.out", n);
  Kokkos::deep_copy(dx, detail::host_ptr_view<const C>(x, n));
  Kokkos::deep_copy(dy, detail::host_ptr_view<const C>(y, n));
  encode2<Bits>(dx, dy, dout);
  Kokkos::deep_copy(detail::host_ptr_view<K>(out, n), dout);
}

template <unsigned Bits>
void encode3_host(const typename Morton<3, Bits>::coord_type* x,
                  const typename Morton<3, Bits>::coord_type* y,
                  const typename Morton<3, Bits>::coord_type* z,
                  typename Morton<3, Bits>::code_type* out, std::size_t n) {
  using C = typename Morton<3, Bits>::coord_type;
  using K = typename Morton<3, Bits>::code_type;
  Kokkos::View<C*> dx("morton.x", n), dy("morton.y", n), dz("morton.z", n);
  Kokkos::View<K*> dout("morton.out", n);
  Kokkos::deep_copy(dx, detail::host_ptr_view<const C>(x, n));
  Kokkos::deep_copy(dy, detail::host_ptr_view<const C>(y, n));
  Kokkos::deep_copy(dz, detail::host_ptr_view<const C>(z, n));
  encode3<Bits>(dx, dy, dz, dout);
  Kokkos::deep_copy(detail::host_ptr_view<K>(out, n), dout);
}

template <unsigned Bits>
void decode2_host(const typename Morton<2, Bits>::code_type* code,
                  typename Morton<2, Bits>::coord_type* x, typename Morton<2, Bits>::coord_type* y,
                  std::size_t n) {
  using C = typename Morton<2, Bits>::coord_type;
  using K = typename Morton<2, Bits>::code_type;
  Kokkos::View<K*> dcode("morton.code", n);
  Kokkos::View<C*> dx("morton.x", n), dy("morton.y", n);
  Kokkos::deep_copy(dcode, detail::host_ptr_view<const K>(code, n));
  decode2<Bits>(dcode, dx, dy);
  Kokkos::deep_copy(detail::host_ptr_view<C>(x, n), dx);
  Kokkos::deep_copy(detail::host_ptr_view<C>(y, n), dy);
}

template <unsigned Bits>
void decode3_host(const typename Morton<3, Bits>::code_type* code,
                  typename Morton<3, Bits>::coord_type* x, typename Morton<3, Bits>::coord_type* y,
                  typename Morton<3, Bits>::coord_type* z, std::size_t n) {
  using C = typename Morton<3, Bits>::coord_type;
  using K = typename Morton<3, Bits>::code_type;
  Kokkos::View<K*> dcode("morton.code", n);
  Kokkos::View<C*> dx("morton.x", n), dy("morton.y", n), dz("morton.z", n);
  Kokkos::deep_copy(dcode, detail::host_ptr_view<const K>(code, n));
  decode3<Bits>(dcode, dx, dy, dz);
  Kokkos::deep_copy(detail::host_ptr_view<C>(x, n), dx);
  Kokkos::deep_copy(detail::host_ptr_view<C>(y, n), dy);
  Kokkos::deep_copy(detail::host_ptr_view<C>(z, n), dz);
}

template <unsigned Dim, unsigned Bits>
void add_host(const typename Morton<Dim, Bits>::code_type* in,
              typename Morton<Dim, Bits>::code_type* out, std::size_t n, unsigned axis,
              long long delta) {
  using K = typename Morton<Dim, Bits>::code_type;
  Kokkos::View<K*> din("morton.in", n), dout("morton.out", n);
  Kokkos::deep_copy(din, detail::host_ptr_view<const K>(in, n));
  add<Dim, Bits>(din, dout, axis, delta);
  Kokkos::deep_copy(detail::host_ptr_view<K>(out, n), dout);
}

}  // namespace kokkos
}  // namespace morton

#endif  // MORTON_KOKKOS_HPP
