#pragma once

#include <cassert>
#include <string>

#include <htslib/hts.h>
#include <htslib/sam.h>


inline std::string_view genomic_substr
(hts_pos_t region_gstart, hts_pos_t from_gpos,
 size_t nchar, std::string_view s)
{
  assert (region_gstart <= from_gpos);
  const auto pos = from_gpos - region_gstart;
  assert (pos < s.size());
  assert ((pos + nchar) < s.size());
  return s.substr (
    static_cast<size_t> (from_gpos - region_gstart),
    nchar
  );
}

// TODO?
std::string to_json (const bam_pileup1_t* p1);

