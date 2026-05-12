#include "pileup-ops.hpp"

#include <filesystem>

namespace pileops {

SamFile open(std::filesystem::path fp) {
  SamFile sf;
  sf.handle.reset (hts_open (fp.c_str(), "w"));

  // NOTE attempting to write a bam1_t
  // to file without having set SQ lines
  // is a segfault if any RNAME is set in the bam1_t
  sf.hdr.reset (sam_hdr_init());
  // BUG placeholder fixed header
  sam_hdr_add_line (sf.hdr.get(), "SQ",
                    "SN", "chr1",
                    "LN", "248956422",  // BUG wrong LN - if pileup data contained pileup params,
                                        // that struct could be used to set this field
                    NULL);
  sam_hdr_write (sf.handle.get(), sf.hdr.get());

  return sf;
}

// header should be extended by this method.
// PileupData should hold the necessary params to do so
// or should also need to pass PileupParams.
int write_pileup(SamFile& sf, const PileupData& pd) {
  const auto read_arr = pd.b1arr.get();
  for (size_t i = 0; i < pd.nread; ++i) {
    const auto rc = sam_write1 (sf.handle.get(), sf.hdr.get(), read_arr + i);
    if (rc < 0) return rc;
  }
  return 0;
}

int write_pileup_auto(const PileupData& pd, std::filesystem::path fp) {
  auto sf = open (fp);
  return write_pileup (sf, pd);
}

}  // end namespace
