/// @file iterate.hpp
/// @brief Region traversal built on the arithmetic core (row-major and Z-order sweeps).
///
/// Two strategies:
///
/// - `for_each_in_box` — row-major (lexicographic) order, advancing with O(1) axis arithmetic so no
///   cell is ever re-encoded from scratch. The cheapest way to sweep a dense axis-aligned region.
/// - `for_each_in_box_zorder` — visits exactly the cells of the region but in Z-order, using the
///   Tropf-Herzog BIGMIN algorithm to skip the gaps between rows. Useful when you need results
///   already sorted by Morton code (range queries on a Z-ordered index).
///
/// SPDX-License-Identifier: MIT

#ifndef MORTON_ITERATE_HPP
#define MORTON_ITERATE_HPP

#include <array>

#include "morton/morton.hpp"

namespace morton {

/// Visit every cell of the inclusive box [lo, hi] in row-major order
/// (axis 0 varies fastest). Advancing from one cell to the next costs O(1)
/// amortised Morton arithmetic; the full coordinate tuple is never re-encoded.
template <unsigned Dim, unsigned Bits, typename Fn>
void for_each_in_box(const std::array<typename Morton<Dim, Bits>::coord_type, Dim>& lo,
                     const std::array<typename Morton<Dim, Bits>::coord_type, Dim>& hi, Fn&& fn) {
  using M = Morton<Dim, Bits>;
  using coord = typename M::coord_type;

  for (unsigned d = 0; d < Dim; ++d)
    if (hi[d] < lo[d])
      return;  // empty box

  M cur = M::encode(lo);
  std::array<coord, Dim> idx = lo;  // track coords without decoding

  for (;;) {
    fn(cur);
    unsigned d = 0;
    for (; d < Dim; ++d) {
      if (idx[d] < hi[d]) {
        ++idx[d];
        cur.inc(d);
        break;
      }
      // axis d exhausted: reset to lo[d] and carry into the next axis
      cur.set(d, lo[d]);
      idx[d] = lo[d];
    }
    if (d == Dim)
      break;  // every axis wrapped -> finished
  }
}

namespace detail {

// Tropf-Herzog BIGMIN: given the Z-order range corners [zmin, zmax] (raw codes
// of the box's lower/upper corners) and a probe code `z` that lies outside the
// box, return the smallest code strictly greater than `z` that is inside the
// box. See H. Tropf & H. Herzog, "Multidimensional Range Search in
// Dynamically Balanced Trees" (1981).
template <unsigned Dim, unsigned Bits>
MORTON_HD typename Morton<Dim, Bits>::code_type bigmin(typename Morton<Dim, Bits>::code_type zmin,
                                                       typename Morton<Dim, Bits>::code_type zmax,
                                                       typename Morton<Dim, Bits>::code_type z) {
  using M = Morton<Dim, Bits>;
  using code = typename M::code_type;

  code result = 0;
  for (int i = int(M::code_bits) - 1; i >= 0; --i) {
    const unsigned d = unsigned(i) % Dim;
    const code bit = code(1) << i;
    const code same_dim_below = M::axis_mask(d) & (bit - 1);

    const int zb = int((z >> i) & 1);
    const int lb = int((zmin >> i) & 1);
    const int hb = int((zmax >> i) & 1);

    if (!zb && !lb && !hb) {
      continue;
    } else if (!zb && !lb && hb) {
      result = (zmin & ~same_dim_below) | bit;  // LOAD(zmin, i, 1)
      zmax = (zmax & ~bit) | same_dim_below;    // LOAD(zmax, i, 0)
    } else if (!zb && lb && hb) {
      return zmin;
    } else if (zb && !lb && !hb) {
      return result;
    } else if (zb && !lb && hb) {
      zmin = (zmin & ~same_dim_below) | bit;  // LOAD(zmin, i, 1)
    } else {                                  // (1,1,1)
      continue;
    }
  }
  return result;
}

// Tropf-Herzog combined LITMAX/BIGMIN. Given Z-order box corners
// [zmin, zmax] and a probe `z`, computes:
//   bigmin = smallest in-box code strictly greater than z, and
//   litmax = largest  in-box code strictly less than z.
// (When the respective value does not exist the result is clamped to a box
// corner; callers that have already established z is outside the box on the
// relevant side get the meaningful answer.)
template <unsigned Dim, unsigned Bits>
MORTON_HD void litmax_bigmin(typename Morton<Dim, Bits>::code_type zmin,
                             typename Morton<Dim, Bits>::code_type zmax,
                             typename Morton<Dim, Bits>::code_type z,
                             typename Morton<Dim, Bits>::code_type& litmax,
                             typename Morton<Dim, Bits>::code_type& bigmin) {
  using M = Morton<Dim, Bits>;
  using code = typename M::code_type;

  litmax = zmin;
  bigmin = zmax;
  for (int i = int(M::code_bits) - 1; i >= 0; --i) {
    const unsigned d = unsigned(i) % Dim;
    const code bit = code(1) << i;
    const code below = M::axis_mask(d) & (bit - 1);  // same-dim bits below i

    const int zb = int((z >> i) & 1);
    const int lb = int((zmin >> i) & 1);
    const int hb = int((zmax >> i) & 1);
    const int v = (zb << 2) | (lb << 1) | hb;

    switch (v) {
      case 0b000:
      case 0b111:
        break;
      case 0b001:
        bigmin = (zmin & ~below) | bit;  // LOAD(zmin, i, 1)
        zmax = (zmax & ~bit) | below;    // LOAD(zmax, i, 0)
        break;
      case 0b011:
        bigmin = zmin;
        return;
      case 0b100:
        litmax = zmax;
        return;
      case 0b101:
        litmax = (zmax & ~bit) | below;  // LOAD(zmax, i, 0)
        zmin = (zmin & ~below) | bit;    // LOAD(zmin, i, 1)
        break;
      default:  // 010, 110 are impossible for zmin <= zmax
        break;
    }
  }
}

}  // namespace detail

/// Smallest Morton code in the inclusive box [lo, hi] that is strictly greater
/// than `probe`. Returns false if there is none.
template <unsigned Dim, unsigned Bits>
bool bigmin_in_box(const std::array<typename Morton<Dim, Bits>::coord_type, Dim>& lo,
                   const std::array<typename Morton<Dim, Bits>::coord_type, Dim>& hi,
                   Morton<Dim, Bits> probe, Morton<Dim, Bits>& out) {
  using M = Morton<Dim, Bits>;
  const auto zmin = M::encode(lo).code();
  const auto zmax = M::encode(hi).code();
  if (probe.code() >= zmax)
    return false;
  typename M::code_type lm, bm;
  detail::litmax_bigmin<Dim, Bits>(zmin, zmax, probe.code(), lm, bm);
  out = M::from_code(bm);
  return true;
}

/// Largest Morton code in the inclusive box [lo, hi] that is strictly less than
/// `probe`. Returns false if there is none.
template <unsigned Dim, unsigned Bits>
bool litmax_in_box(const std::array<typename Morton<Dim, Bits>::coord_type, Dim>& lo,
                   const std::array<typename Morton<Dim, Bits>::coord_type, Dim>& hi,
                   Morton<Dim, Bits> probe, Morton<Dim, Bits>& out) {
  using M = Morton<Dim, Bits>;
  const auto zmin = M::encode(lo).code();
  const auto zmax = M::encode(hi).code();
  if (probe.code() <= zmin)
    return false;
  typename M::code_type lm, bm;
  detail::litmax_bigmin<Dim, Bits>(zmin, zmax, probe.code(), lm, bm);
  out = M::from_code(lm);
  return true;
}

/// Test whether `m` lies inside the inclusive box [lo, hi].
template <unsigned Dim, unsigned Bits>
MORTON_HD bool in_box(const std::array<typename Morton<Dim, Bits>::coord_type, Dim>& lo,
                      const std::array<typename Morton<Dim, Bits>::coord_type, Dim>& hi,
                      Morton<Dim, Bits> m) {
  for (unsigned d = 0; d < Dim; ++d) {
    auto c = m.get(d);
    if (c < lo[d] || c > hi[d])
      return false;
  }
  return true;
}

/// Visit every cell of the inclusive box [lo, hi] in Z-order (increasing
/// Morton code), skipping the gaps between rows with BIGMIN.
template <unsigned Dim, unsigned Bits, typename Fn>
void for_each_in_box_zorder(const std::array<typename Morton<Dim, Bits>::coord_type, Dim>& lo,
                            const std::array<typename Morton<Dim, Bits>::coord_type, Dim>& hi,
                            Fn&& fn) {
  using M = Morton<Dim, Bits>;
  using code = typename M::code_type;

  for (unsigned d = 0; d < Dim; ++d)
    if (hi[d] < lo[d])
      return;

  const code zmin = M::encode(lo).code();
  const code zmax = M::encode(hi).code();

  code z = zmin;
  for (;;) {
    M m = M::from_code(z);
    if (in_box<Dim, Bits>(lo, hi, m)) {
      fn(m);
      if (z == zmax)
        break;
      ++z;
    } else {
      code lm, next;
      detail::litmax_bigmin<Dim, Bits>(zmin, zmax, z, lm, next);
      if (next <= z)
        break;  // safety: guarantee forward progress
      z = next;
    }
  }
}

}  // namespace morton

#endif  // MORTON_ITERATE_HPP
