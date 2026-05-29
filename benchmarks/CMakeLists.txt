add_executable(morton_bench bench.cpp)
target_link_libraries(morton_bench PRIVATE morton)
target_include_directories(morton_bench PRIVATE ${CMAKE_SOURCE_DIR}/third_party)
target_compile_options(morton_bench PRIVATE -O3)
