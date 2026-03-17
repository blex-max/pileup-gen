#pragma once

#include <cstddef>
#include <functional>
#include <vector>

#include <htslib/sam.h>

#include "read-ops.hpp"

// I'm not yet clear enough on what
// good decisions and bad decisions
// would be on the api for generating
// reads from one or more "fuzzy" specs.
// So I'm continuting with exact targets for now.

// TODO obvious patterns to implement:
// n target exactly specified
// n target with some fields exactly specified and some random (or sequential)
// n target and n random
// n.b must support paired end
// (NOTE later, this can grow into an additional API not for specific targets but for
// more general model based simulation like "poorly aligning reads" and so on)

using ReadV = std::vector<bam1_t*>;
using ReadArgV = std::vector<std::pair<size_t, readops::ReadData>>;


ReadV generate_reads (const ReadArgV& rav);


ReadV gen_fuzzy_read (std::string_view ref, size_t read_len);
