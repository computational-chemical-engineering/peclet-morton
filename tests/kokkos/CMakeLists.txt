# Portable Kokkos backend test. Built only when MORTON_ENABLE_KOKKOS=ON (the
# parent gates this add_subdirectory on it). Links Kokkos::kokkos + the define
# via the morton interface target. Device sources compile as CXX through Kokkos
# (Kokkos 5.x wires the device flags onto any target linking Kokkos::kokkos), so
# this stays a plain .cpp -- see ../../cmake/SuiteKokkos.cmake.
add_executable(morton_kokkos_tests test_kokkos.cpp)
target_link_libraries(morton_kokkos_tests PRIVATE morton)

add_test(NAME morton_kokkos_tests COMMAND morton_kokkos_tests)
