# vcpkg port for the header-only morton-arithmetic library.
#
# Before submitting to the vcpkg registry, set REF to the release tag and fill
# SHA512 with the value reported by `vcpkg install morton` on first attempt
# (or `vcpkg_from_github(... SHA512 0)` once and copy the expected hash).

vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO computational-chemical-engineering/morton_artithmetic
    REF v0.1.0
    SHA512 0
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
