#include "pileup-ops.hpp"

#include <filesystem>

namespace pileops {

int write_pileup (const PileupData& pd, std::filesystem::path fp) {
  // BUG doesn't make sense to be returning different error ints...

  auto sam_handle = hts_open (fp.c_str(), "w");

  // BUG placeholder fixed header
  auto hdr = sam_hdr_init ();
  // NOTE attempting to write a bam1_t
  // to file without having set SQ lines
  // is a segfault if any RNAME
  // is set in the bam1_t
  const std::string tid ("chr1");
  sam_hdr_add_line (hdr, "SQ",
                   "SN", tid.c_str(),
                   "LN", "248956422",  // BUG wrong ln - if pileup data contained pileup params could use that
                   NULL);
  if (const auto rc = sam_hdr_write (sam_handle, hdr);
      rc < 0)
  {
        return rc;
  };

  const auto read_arr = pd.b1arr.get();
  for (size_t i = 0; i < pd.nread; ++i) {
    const auto rc = sam_write1(sam_handle, hdr, read_arr + i);
    if (rc < 0) {
      return rc;
    }
  }

  // I really need unique ptrs
  hts_flush(sam_handle);
  hts_close(sam_handle);
  sam_hdr_destroy(hdr);

  return 0;
}

}  // end namespace
