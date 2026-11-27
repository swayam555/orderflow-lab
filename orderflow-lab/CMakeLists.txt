add_library(ofl STATIC
  src/order_book.cpp src/reference_book.cpp)
target_include_directories(ofl PUBLIC include)
target_link_libraries(ofl PUBLIC qc)
qc_target_defaults(ofl)

qc_add_test(ofl.test_book_basic tests/test_book_basic.cpp "unit")
qc_add_test(ofl.test_book_random tests/test_book_random.cpp "unit")
qc_add_test(ofl.test_book_random_slow tests/test_book_random_slow.cpp "slow")
target_include_directories(ofl_test_book_random_slow PRIVATE tests)

add_executable(ofl_bench_engine bench/ofl_bench_engine.cpp)
qc_target_defaults(ofl_bench_engine)
target_link_libraries(ofl_bench_engine PRIVATE ofl qc)
target_compile_definitions(ofl_bench_engine PRIVATE OFL_RESULTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/results")
add_test(NAME ofl.bench_engine_verify COMMAND ofl_bench_engine --ops=300000 --warmup=100000 --reps=1 --verify=1
         --out=${CMAKE_BINARY_DIR}/bench_tmp WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(ofl.bench_engine_verify PROPERTIES LABELS "bench")
