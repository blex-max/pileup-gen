#include <filesystem>
#include <memory>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "generate-pileup.hpp"

namespace pileops {

struct HtsFileDeleter {
  void operator()(htsFile* f) const { hts_flush (f); hts_close (f); }
};

struct SamHdrDeleter {
  void operator()(sam_hdr_t* h) const { sam_hdr_destroy (h); }
};

// hdr declared before handle: destruction is handle first (flush+close), then hdr
struct SamFile {
  std::unique_ptr<sam_hdr_t, SamHdrDeleter> hdr;
  std::unique_ptr<htsFile, HtsFileDeleter> handle;
};

SamFile open(std::filesystem::path fp);
int write_pileup(SamFile& sf, const PileupData& pd);

int write_pileup_auto(const PileupData& pd, std::filesystem::path fp);

}
