// tests-dynamic_bitset.cpp
#include "catch.hpp"

#include "libnormaliz/dynamic_bitset.h"

TEST_CASE("resize, size, empty", "[dynamic_bitset]") {

    libnormaliz::dynamic_bitset b(5);

    REQUIRE(b.size() == 5);
    REQUIRE_FALSE(b.empty());
    REQUIRE(b.none());

    SECTION("resizing bigger changes size") {
        b.resize(10);

        REQUIRE(b.size() == 10);
		REQUIRE_FALSE(b.empty());
        REQUIRE(b.none());
    }
    SECTION("resizing smaller changes size") {
        b.resize(0);

        REQUIRE(b.size() == 0);
		REQUIRE(b.empty());
        REQUIRE(b.none());
    }
}

TEST_CASE("find_first, find_next", "[dynamic_bitset]") {

	SECTION("empty-bitset") {
		libnormaliz::dynamic_bitset b;

		REQUIRE(b.find_first() == b.npos);
		REQUIRE(b.find_next(0) == b.npos);
		REQUIRE(b.find_next(1) == b.npos);
		REQUIRE(b.find_next(100) == b.npos);
		REQUIRE(b.find_next(b.npos) == b.npos);
	}
 
	SECTION("non-empty bitset") {
	    // TODO: use Catch2 RandomIntegerGenerator to generate a bunch of
	    // random bitsets and test with them
		libnormaliz::dynamic_bitset b(101);

		REQUIRE(b.find_first() == b.npos);
		REQUIRE(b.find_next(0) == b.npos);
		REQUIRE(b.find_next(1) == b.npos);
		REQUIRE(b.find_next(9) == b.npos);
		REQUIRE(b.find_next(10) == b.npos);
		REQUIRE(b.find_next(11) == b.npos);
		REQUIRE(b.find_next(99) == b.npos);
		REQUIRE(b.find_next(100) == b.npos);
		REQUIRE(b.find_next(101) == b.npos);
		REQUIRE(b.find_next(b.npos) == b.npos);

		b.set(100);
		REQUIRE(b.find_first() == 100);
		REQUIRE(b.find_next(0) == 100);
		REQUIRE(b.find_next(1) == 100);
		REQUIRE(b.find_next(9) == 100);
		REQUIRE(b.find_next(10) == 100);
		REQUIRE(b.find_next(11) == 100);
		REQUIRE(b.find_next(99) == 100);
		REQUIRE(b.find_next(100) == b.npos);
		REQUIRE(b.find_next(101) == b.npos);
		REQUIRE(b.find_next(b.npos) == b.npos);

		b.set(10);
		REQUIRE(b.find_first() == 10);
		REQUIRE(b.find_next(0) == 10);
		REQUIRE(b.find_next(1) == 10);
		REQUIRE(b.find_next(9) == 10);
		REQUIRE(b.find_next(10) == 100);
		REQUIRE(b.find_next(11) == 100);
		REQUIRE(b.find_next(99) == 100);
		REQUIRE(b.find_next(100) == b.npos);
		REQUIRE(b.find_next(101) == b.npos);
		REQUIRE(b.find_next(b.npos) == b.npos);

		b.set(1);
		REQUIRE(b.find_first() == 1);
		REQUIRE(b.find_next(0) == 1);
		REQUIRE(b.find_next(1) == 10);
		REQUIRE(b.find_next(9) == 10);
		REQUIRE(b.find_next(10) == 100);
		REQUIRE(b.find_next(11) == 100);
		REQUIRE(b.find_next(99) == 100);
		REQUIRE(b.find_next(100) == b.npos);
		REQUIRE(b.find_next(101) == b.npos);
		REQUIRE(b.find_next(b.npos) == b.npos);
	}

}


TEST_CASE("all, any, count", "[dynamic_bitset]") {

	SECTION("empty-bitset") {
		libnormaliz::dynamic_bitset b;

		REQUIRE_FALSE(b.any());
		REQUIRE(b.none());
		REQUIRE(b.count() == 0);
	}

	SECTION("non-empty bitset") {
	    // TODO: use Catch2 RandomIntegerGenerator to generate a bunch of
	    // random bitsets and test with them
		libnormaliz::dynamic_bitset b(101);

		REQUIRE_FALSE(b.any());
		REQUIRE(b.none());
		REQUIRE(b.count() == 0);

		b.set(100);
		REQUIRE(b.any());
		REQUIRE_FALSE(b.none());
		REQUIRE(b.count() == 1);

		b.set(10);
		REQUIRE(b.any());
		REQUIRE_FALSE(b.none());
		REQUIRE(b.count() == 2);

		b.set(1);
		REQUIRE(b.any());
		REQUIRE_FALSE(b.none());
		REQUIRE(b.count() == 3);
	}

}
