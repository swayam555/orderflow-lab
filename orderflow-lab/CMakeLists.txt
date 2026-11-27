add_library(ofl STATIC
  src/order_book.cpp)
target_include_directories(ofl PUBLIC include)
target_link_libraries(ofl PUBLIC qc)
qc_target_defaults(ofl)

qc_add_test(ofl.test_book_basic tests/test_book_basic.cpp "unit")
