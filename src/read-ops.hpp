#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <variant>
#include <vector>

#include <htslib/sam.h>


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

using AuxTag = const char[3];
using AuxData = std::variant<int64_t, float, std::string>;
int append_aux (bam1_t* b, AuxTag name, AuxData data);
// using AuxArrayData... todo


struct ReadSpec {
  std::string qseq{};   // emtpy == NULL
  std::string qqual{};  // empty == NULL
  std::string qname{};  // empty == NULL
  CigV qcig{};      // empty == unaligned
  hts_pos_t lmost_pos=0;
  hts_pos_t mate_lmost_pos=0;
  uint16_t flag=0;
  int32_t tid=0;
  int32_t mate_tid=-1;
  uint8_t mapq=0;
  // TLEN/isize left out at least for now
  // as the field is problematic and nonstandard
  // (per discussion with sam team).
  std::map<AuxTag, AuxData> aux;
  // auxarray...
};
std::ostream& operator<< (std::ostream& os, const ReadSpec& rs);  // serialise


// wrapper for bam_set1
int set_bam1 (const ReadSpec& rs, bam1_t* b);


// of undetermined utility
// void ins (bam1_t* b, size_t pos, std::string ins);
// void del (bam1_t* b, size_t pos, size_t len);
// void sub (bam1_t* b, size_t pos, char base);
// void clip (bam1_t* b);

// void assign_qual (bam1_t* b, QualModel qual);

}  // end namespace
