add_executable(morton_tests
    test_main.cpp
    test_encode_decode.cpp
    test_arithmetic.cpp
    test_iterate.cpp
    test_constexpr.cpp
    test_wide.cpp
    test_neighbors.cpp
    test_batch.cpp)
target_link_libraries(morton_tests PRIVATE morton)
target_include_directories(morton_tests PRIVATE ${CMAKE_SOURCE_DIR}/third_party)

add_test(NAME morton_tests COMMAND morton_tests)
