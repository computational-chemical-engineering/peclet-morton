# vcpkg port for the header-only peclet-morton library.
#
# Before submitting to the vcpkg registry, set REF to the release tag and fill
# SHA512 with the value reported by `vcpkg install morton` on first attempt
# (or `vcpkg_from_github(... SHA512 0)` once and copy the expected hash).

vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO computational-chemical-engineering/peclet-morton
    REF v1.0.2
    SHA512 c2a67f09e88d56bec760653f33ea896988cce0be319e4c0dd0bf304ffc211d3ae1d07a9867aac8108d18d05855e68fe65c3f607a5361131e8329d916741cd839
    HEAD_REF main
)

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        -DMORTON_BUILD_TESTS=OFF
        -DMORTON_BUILD_BENCHMARKS=OFF
        -DMORTON_BUILD_BINDINGS=OFF
)

vcpkg_cmake_install()

vcpkg_cmake_config_fixup(PACKAGE_NAME morton CONFIG_PATH lib/cmake/morton)

# Header-only: drop the empty debug tree.
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug")

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE")
