#include <cassert>
#include <htslib/sam.h>

#include "generate-reads.hpp"

ReadV generate_reads (const ReadArgV& rav) {
    ReadV out;
    for (const auto& [n, args] : rav) {
        for (size_t i = 0; i < n; ++i) {
            out.push_back (readops::create_read (args));
        }
    }
    return out;
}

