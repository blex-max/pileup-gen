#include <string>

#include <htslib/hts.h>
#include <htslib/sam.h>


inline std::string get_seq (const bam1_t *b) {
  std::string out;

  auto n = static_cast<size_t> (b->core.l_qseq);
  auto seq8 = bam_get_seq(b);

  for (size_t i=0; i < n; ++i) {
    out.push_back(seq_nt16_str[bam_seqi(seq8, i)]);
  }

  return out;
}
