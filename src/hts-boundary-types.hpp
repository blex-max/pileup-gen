#pragma once

#include <htslib/hts.h>
#include <htslib/sam.h>


namespace hts {

struct SamOut {
  samFile* ofp;
  sam_hdr_t* hdr;
};



}  // end namespace
