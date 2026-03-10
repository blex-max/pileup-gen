#pragma once

#include <format>
#include <iostream>
#include <stdexcept>
#include <string>

#include <htslib/sam.h>
#include <vector>


namespace readops {

// NOTE
// create this functionality first,
// then streamline once working
// - there's obviously some overlap between
// structs/common properties between functions
using CigarVec = std::vector<std::pair<size_t, char>>;
struct ReadArgs {
  std::string qseq{};   // emtpy == NULL
  std::string qqual{};  // empty == NULL
  std::string qname{};  // empty == NULL
  CigarVec qcig{};      // empty == unaligned
  hts_pos_t lmost_pos=0;
  hts_pos_t mate_lmost_pos=0;
  uint16_t flag=0;
  int32_t tid=0;
  int32_t mate_tid=0;
  uint8_t mapq=0;
  // TLEN/isize left out at least for now
  // as the field is problematic and nonstandard
  // (per discussion with sam team).
  // l_aux left out *for now*, under further consideration
};

// wrapper for bam_set1
inline bam1_t* create_read (const ReadArgs& ra) {
  auto b = bam_init1();

  std::string cig_s;
  size_t cig_nop = ra.qcig.size();
  if (cig_nop > 0) {
    cig_s.reserve(cig_nop);
    for (const auto cp : ra.qcig) {
      cig_s.append(std::format ("{}{}", cp.first, cp.second));
    }
    // TODO check number processed (the return of *_parse_cigar)
  } else {
    cig_s = "*";
  }
  uint32_t  *cig_buf = static_cast<uint32_t *> (
      malloc (cig_nop * sizeof (uint32_t))
  );
  size_t buf_alloc;
  sam_parse_cigar (cig_s.c_str(), NULL, &cig_buf, &buf_alloc);

  const auto rc = bam_set1 (
      b,
      ra.qname.size(),
      ra.qname.empty() ? NULL : ra.qname.c_str(),
      cig_nop > 0 ? ra.flag : ra.flag | BAM_FUNMAP ,
      -1,
      ra.lmost_pos,
      ra.mapq,
      cig_nop,
      cig_nop > 0 ? cig_buf : NULL,  // cigar set later
      -1,
      ra.mate_lmost_pos,
      0,  // isize ignored for now
      ra.qseq.size(),
      ra.qseq.empty() ? NULL : ra.qseq.c_str(),
      ra.qqual.empty() ? NULL : ra.qqual.c_str(),
      0  // TODO, perhaps
  );

  if (rc < 0) {
    throw std::runtime_error ("failed to set bam1_t");
  }

  return b;
};


// probably unnecessary
void ins (bam1_t* b, size_t pos, std::string ins);
void del (bam1_t* b, size_t pos, size_t len);
void sub (bam1_t* b, size_t pos, char base);
// void clip (bam1_t* b);

// void assign_qual (bam1_t* b, QualModel qual);

}
