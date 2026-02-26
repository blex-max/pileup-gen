#pragma once

#include <cstddef>
#include <random>
#include <vector>

#include <htslib/sam.h>

// which fields to generate
// or keep fixed
// std::variant? std::optional?
// sub-struct?
// in any case sort of like a generic
// template
struct ReadProps {
  
};


using ReadVec = std::vector<bam1_t*>;
using ReadSpecVec = std::vector<std::pair<ReadProps, size_t>>;
class ReadGen {
  private:
  std::mt19937 rng();

  public:
  ReadProps rp;

  ReadVec
  generate_reads (ReadSpecVec rsv);
};
