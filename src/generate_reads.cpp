#include "generate-reads.hpp"
#include <htslib/sam.h>

// So far so good
// but there's an even simpler grain
// to figure out first
// generating an exactly specified target read - TODO
ReadVec ReadGen::generate_reads (size_t n) {
  ReadVec out(n);  // default construct bam1_t objects

  // create partially set bam1_t with templated properties
  // as set by read props
  // (TODO, placeholder below)
  std::string qseq = "AAAA";
  const auto ncigop  = 1;
  std::string cs{"*"};
  auto b_tmpl = bam_init1();
  uint32_t  *cig_buf = static_cast<uint32_t *> (
      malloc (ncigop * sizeof (uint32_t))
  );
  size_t buf_alloc;
  sam_parse_cigar (cs.c_str(), NULL, &cig_buf, &buf_alloc);
  // bam_parse_cigar (cs.c_str(), NULL, b);
  bam_set1 (
      b_tmpl,
      0,
      NULL,
      0,
      0,
      0,
      0,
      ncigop,
      cig_buf,
      0,
      0,
      0,
      qseq.size(),
      qseq.c_str(),
      NULL,
      0
  );

  for (auto& b : out) {
    // bam_set1 (
    //   &b
    //   b_tmpl.l_qname...
    //   next_qname()
    //   and so on
    // )
  }

  return out;
}
