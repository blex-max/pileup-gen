#include "read-ops.hpp"

#include <cstdint>
#include <format>
#include <ostream>
#include <ios>

#include <htslib/sam.h>


namespace readops {

// ReadSpec
std::ostream& operator<< (std::ostream& os, const ReadSpec& rs) {
  const auto s = std::format(
    "qname={}\nqseq={}\nqqual={}\nqcig={}\n"
    "flag={}\ntid={}\nstart={}\nmapq={}\nmate_tid={}\nmate_start={}",
    rs.qname, rs.qseq, "STRINGIFY-NOT-IMPLEMETED", "STRINGIFY-NOT-IMPLEMENTED",
    rs.flag, rs.tid, rs.lmost_pos, rs.mapq, rs.mate_tid, rs.mate_lmost_pos
  );
  os.write(s.data(), static_cast<std::streamsize>(s.size()));
  return os;
}


// functor helper for std::visit
template<class... Ts>
struct overloads : Ts... { using Ts::operator()...; };
int append_aux
(bam1_t* b, AuxTag name, AuxData data)
{
  if (name.size() != 2) return -1;
  const auto vistor = overloads {
    [b, name] (std::string s) {
      return bam_aux_update_str (b, name.c_str(), -1, s.c_str());
    },
    [b, name] (int64_t i) {
      return bam_aux_update_int (b, name.c_str(), i);
    },
    [b, name] (float f) {
      return bam_aux_update_float (b, name.c_str(), f);
    }
  };

  return std::visit (vistor, data);
}


int set_bam1 (const ReadSpec& rs, bam1_t* b) {

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
      rs.tid,
      rs.lmost_pos,
      rs.mapq,
      cig_nop,
      cig_nop > 0 ? cig_arr : NULL,  // cigar set later
      rs.mate_tid,
      rs.mate_lmost_pos,
      0,  // isize ignored for now
      rs.qseq.size(),
      rs.qseq.empty() ? NULL : rs.qseq.c_str(),
      rs.qqual.empty() ? NULL : rs.qqual.c_str(),
      0  // append later
  );

  delete[] cig_arr;

  if (rc < 0) {
    return rc;
  }

  for (const auto& [name, data] : rs.aux) {
    append_aux (b, name, data);
  }

  return 0;
};


}  // end namespace readops
