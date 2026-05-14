#include "accessors.hpp"

#include <htslib/hts.h>
#include <string>

#include <htslib/sam.h>

namespace htsacc {


hts_pos_t qpos (const bam_pileup1_t* p1) {
  return p1->qpos;
}


hts_pos_t gstart (const bam_pileup1_t* p1) {
  return p1->b->core.pos;
}

char pileup_base (const bam_pileup1_t* p1) {
  return seq_nt16_str[bam_seqi(bam_get_seq(p1->b), p1->qpos)];
}

uint16_t flag (const bam_pileup1_t* p1) {
  return p1->b->core.flag;
}


std::string seq (const bam1_t *b, size_t qpos, size_t n)
{
  if (n == 0 || (qpos + n) > b->core.l_qseq) {
    n = b->core.l_qseq;
  }
  std::string out;
  out.resize(n);
  auto seq_nib = bam_get_seq(b);

  for (size_t i=0; i < n; ++i) {
    out.push_back(seq_nt16_str[bam_seqi(seq_nib, qpos + i)]);
  }

  return out;
}
std::string seq (const bam_pileup1_t* p1, size_t qpos, size_t n)
{
  return seq (p1->b, qpos, n);
}



std::string get_seq_genomic (const bam_pileup1_t* p1, size_t gpos, size_t n=0) {
  const auto qpos = gpos - p1->b->core.pos;
  return seq (p1, qpos, n);
}

uint8_t base_qual (const bam_pileup1_t *p1) {
  const auto qpos = p1->qpos;

  const auto bq = bam_get_qual(p1->b);  // get qual string arr
  if (bq == NULL) {
    return 0;
  }
  return *(bq + qpos);
}

}  // end namespace


