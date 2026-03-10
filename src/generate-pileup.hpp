#include <htslib/hts.h>
#include <htslib/sam.h>
#include <random>


enum struct BaseEvents : size_t {
  rand,
  del,
  A,
  C,
  T,
  G,
  // other codes TODO
};


// concrete events
struct EventSpec {
  EventSpec () = delete;
  BaseEvents base;
  int qpos=-1;        // <0 random
  int indel=0;        // <0 del, >0 ins, length indel
  std::string ins{};  // empty is random, mismatched length is err
};
// NOTE the above *could*, in principle, describe any base in a read


// NOTE WIP
using PileupVec = std::vector<bam_pileup1_t>;
class PileupGen {
  private:
  std::mt19937 rng;

  // can provide any single, and the others will be derived
  std::string ref_region;  // ref pos is center, therefore N must be odd
  size_t ref_pos;          // see above
  size_t read_len;         // <= (ref_region.size() + 1) / 2

  public:
  PileupGen () = delete;

  // random
  PileupVec
  generate (double pt_non_ref, size_t nreads);

  // percentage spec
  PileupVec
  generate (const std::vector<std::pair<double, EventSpec>> &events, size_t nreads);

  // explicit specification
  PileupVec
  generate (const std::vector<std::pair<size_t, EventSpec>> &events);

};

// TODO
// generate_pileup.cpp
// using create_read
