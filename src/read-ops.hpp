#pragma once

#include <cstdint>
#include <string>

#include <htslib/sam.h>
#include <vector>


namespace readops {

enum cigarcode : int {
  match = BAM_CMATCH,
  equal = BAM_CEQUAL,
  diff = BAM_CDIFF,
  del = BAM_CDEL,
  ins = BAM_CINS,
  softclip = BAM_CSOFT_CLIP,
  // other codes unused
};

using CigV = std::vector<std::pair<size_t, cigarcode>>;
struct ReadSpec {
  std::string qseq{};   // emtpy == NULL
  std::string qqual{};  // empty == NULL
  std::string qname{};  // empty == NULL
  CigV qcig{};      // empty == unaligned
  hts_pos_t lmost_pos=0;
  hts_pos_t mate_lmost_pos=0;
  uint16_t flag=0;
  int32_t tid=0;
  int32_t mate_tid=0;
  uint8_t mapq=0;
  // TLEN/isize left out at least for now
  // as the field is problematic and nonstandard
  // (per discussion with sam team).
  // aux left out for now, due further consideration
};

// wrapper for bam_set1
inline int set_bam1 (const ReadSpec& rs, bam1_t* b) {

  size_t cig_nop = rs.qcig.size();
  uint32_t* cig_arr = static_cast<uint32_t*> (new uint32_t[cig_nop]);
  if (cig_nop > 0) {
    for (size_t opi = 0; opi < cig_nop; ++opi) {
      const auto op = rs.qcig[opi];
      cig_arr[opi] = bam_cigar_gen (op.first, op.second);
    }
  }

  const auto rc = bam_set1 (
      b,
      rs.qname.size(),
      rs.qname.empty() ? NULL : rs.qname.c_str(),
      cig_nop > 0 ? rs.flag : rs.flag | BAM_FUNMAP ,
      -1,
      rs.lmost_pos,
      rs.mapq,
      cig_nop,
      cig_nop > 0 ? cig_arr : NULL,  // cigar set later
      -1,
      rs.mate_lmost_pos,
      0,  // isize ignored for now
      rs.qseq.size(),
      rs.qseq.empty() ? NULL : rs.qseq.c_str(),
      rs.qqual.empty() ? NULL : rs.qqual.c_str(),
      0  // TODO, perhaps
  );

  delete cig_arr;

  return rc;
};


// of undetermined utility
// void ins (bam1_t* b, size_t pos, std::string ins);
// void del (bam1_t* b, size_t pos, size_t len);
// void sub (bam1_t* b, size_t pos, char base);
// void clip (bam1_t* b);

// void assign_qual (bam1_t* b, QualModel qual);

}
