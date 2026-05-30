add_executable(morton_bench bench.cpp)
target_link_libraries(morton_bench PRIVATE morton)
target_include_directories(morton_bench PRIVATE ${CMAKE_SOURCE_DIR}/third_party)
target_compile_options(morton_bench PRIVATE -O3)

add_executable(morton_bench_batch bench_batch.cpp)
target_link_libraries(morton_bench_batch PRIVATE morton)
target_compile_options(morton_bench_batch PRIVATE -O3)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|X86_64|amd64|AMD64|i[3-6]86)$")
    target_compile_options(morton_bench_batch PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>:-mavx2>)
endif()
