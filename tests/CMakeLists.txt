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

# Contract: the plain BMI2=OFF build (no -mbmi2, no runtime dispatch) emits no PDEP/PEXT at all,
# so it can never SIGILL on a CPU without BMI2. Checked by disassembling the test binary.
if(NOT MORTON_ENABLE_BMI2 AND NOT MORTON_ENABLE_RUNTIME_DISPATCH)
    find_program(MORTON_OBJDUMP NAMES objdump llvm-objdump)
    if(MORTON_OBJDUMP)
        add_test(NAME morton_no_pdep_pext
                 COMMAND ${CMAKE_COMMAND} -DOBJDUMP=${MORTON_OBJDUMP}
                         -DBINARY=$<TARGET_FILE:morton_tests>
                         -P ${CMAKE_CURRENT_SOURCE_DIR}/check_no_pdep.cmake)
    else()
        message(STATUS "morton: objdump not found; the no-pdep/pext contract test is not registered")
    endif()
endif()
