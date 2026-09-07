# ctest helper: fail if the given binary contains any PDEP/PEXT instruction.
# Usage: cmake -DOBJDUMP=<objdump> -DBINARY=<file> -P check_no_pdep.cmake
execute_process(COMMAND "${OBJDUMP}" -d --no-show-raw-insn "${BINARY}"
                OUTPUT_VARIABLE _asm RESULT_VARIABLE _rc ERROR_VARIABLE _err)
if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "objdump failed on ${BINARY}: ${_err}")
endif()
string(REGEX MATCHALL "[ \t](pdep|pext)[ \t]" _hits "${_asm}")
list(LENGTH _hits _n)
if(_n GREATER 0)
    message(FATAL_ERROR "${_n} pdep/pext instruction(s) in ${BINARY}: the MORTON_ENABLE_BMI2=OFF "
                        "build must be free of them (software path only)")
endif()
message(STATUS "no pdep/pext in ${BINARY}")
