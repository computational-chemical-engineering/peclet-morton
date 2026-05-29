#include "doctest.h"
#include "morton/iterate.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>

using namespace morton;

TEST_CASE("2D box iteration (row-major and Z-order) matches brute force") {
    using M = Morton<2, 5>;  // 32x32
    using coord = M::coord_type;

    for (coord lx = 0; lx < 32; lx += 3)
        for (coord ly = 0; ly < 32; ly += 3)
            for (coord hx = lx; hx < 32; hx += 5)
                for (coord hy = ly; hy < 32; hy += 5) {
                    std::array<coord, 2> lo{lx, ly}, hi{hx, hy};

                    std::vector<std::uint32_t> expect;
                    for (coord y = ly; y <= hy; ++y)
                        for (coord x = lx; x <= hx; ++x)
                            expect.push_back(M::encode(x, y).code());
                    std::sort(expect.begin(), expect.end());

                    std::vector<std::uint32_t> rowmajor;
                    for_each_in_box<2, 5>(lo, hi, [&](M m) { rowmajor.push_back(m.code()); });
                    auto sorted = rowmajor;
                    std::sort(sorted.begin(), sorted.end());
                    CHECK(sorted == expect);

                    std::vector<std::uint32_t> zorder;
                    for_each_in_box_zorder<2, 5>(lo, hi, [&](M m) { zorder.push_back(m.code()); });
                    CHECK(std::is_sorted(zorder.begin(), zorder.end()));
                    CHECK(zorder == expect);
                }
}

TEST_CASE("3D Z-order iteration matches brute force") {
    using M = Morton<3, 4>;  // 16^3
    using coord = M::coord_type;
    std::array<coord, 3> lo{2, 3, 1}, hi{9, 7, 6};

    std::vector<std::uint64_t> expect;
    for (coord z = lo[2]; z <= hi[2]; ++z)
        for (coord y = lo[1]; y <= hi[1]; ++y)
            for (coord x = lo[0]; x <= hi[0]; ++x)
                expect.push_back(M::encode(x, y, z).code());
    std::sort(expect.begin(), expect.end());

    std::vector<std::uint64_t> zorder;
    for_each_in_box_zorder<3, 4>(lo, hi, [&](M m) { zorder.push_back(m.code()); });
    CHECK(std::is_sorted(zorder.begin(), zorder.end()));
    CHECK(zorder == expect);
}

TEST_CASE("degenerate boxes") {
    using M = Morton<2, 6>;
    using coord = M::coord_type;

    // single cell
    int count = 0;
    for_each_in_box<2, 6>({5, 7}, {5, 7}, [&](M) { ++count; });
    CHECK(count == 1);

    // empty box (hi < lo) visits nothing
    count = 0;
    for_each_in_box<2, 6>({5, 7}, {4, 7}, [&](M) { ++count; });
    CHECK(count == 0);

    count = 0;
    for_each_in_box_zorder<2, 6>({5, 7}, {4, 7}, [&](M) { ++count; });
    CHECK(count == 0);

    (void)sizeof(coord);
}
