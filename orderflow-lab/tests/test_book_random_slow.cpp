// Same differential test at scale (label slow): 50 sequences x 200000 ops.
#define OFL_SLOW_RANDOM 1
#include "test_book_random.cpp"

TEST_CASE("randomized differential vs reference book (slow: 50 x 200000 ops)") {
    for (std::uint64_t s = 0; s < 50; ++s) {
        const SeqResult r = run_sequence(s, 200000, 1000);
        CHECK_MESSAGE(r.ok, r.msg);
        if (!r.ok) break;
    }
}
