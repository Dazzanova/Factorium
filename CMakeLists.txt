# Factorium — Phase 1: core engine and correctness.
cmake_minimum_required(VERSION 3.16)
project(factorium C)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

add_library(factorium_core STATIC
  src/common.c
  src/modular_arithmetic.c
  src/primality.c
  src/fermat.c
  # src/pollard_rho.c
  # src/pollard_p_minus_one.c
  # src/factorization.c
)
target_include_directories(factorium_core PUBLIC include)

# add_executable(factorium src/main.c)
# target_link_libraries(factorium PRIVATE factorium_core)

# Strict warnings on everything we build.
set(FACTORIUM_WARNING_FLAGS -Wall -Wextra -Wpedantic -Wshadow
                            -Wstrict-prototypes -Wmissing-prototypes)
foreach(t factorium_core)
  target_compile_options(${t} PRIVATE ${FACTORIUM_WARNING_FLAGS})
endforeach()

option(FACTORIUM_SANITIZERS "Build tests with ASan+UBSan (development)" OFF)

enable_testing()

set(FACTORIUM_TESTS
  test_modular_arithmetic
  test_primality
  test_fermat
  # test_pollard_rho
  # test_pollard_p_minus_one
  # test_factorization
  # test_cli
)

foreach(name IN LISTS FACTORIUM_TESTS)
  add_executable(${name} tests/${name}.c)
  target_link_libraries(${name} PRIVATE factorium_core)
  target_compile_options(${name} PRIVATE ${FACTORIUM_WARNING_FLAGS})
  if(FACTORIUM_SANITIZERS)
    target_compile_options(${name} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(${name} PRIVATE -fsanitize=address,undefined)
  endif()
  if(name STREQUAL "test_cli")
    # The CLI test receives the factorium binary path as argv[1].
    add_test(NAME cli COMMAND ${name} $<TARGET_FILE:factorium>)
  else()
    add_test(NAME ${name} COMMAND ${name})
  endif()
endforeach()

# Optional: validate JSON output with Python 3 (stdlib only) when available.
# find_program(PYTHON3 NAMES python3)
# if(PYTHON3)
#   add_test(NAME cli_json
#            COMMAND ${PYTHON3} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_cli_json.py
#                    $<TARGET_FILE:factorium>)
# endif()
