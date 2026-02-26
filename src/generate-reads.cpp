#include "generate-reads.hpp"
#include <cassert>
#include <htslib/sam.h>

using ReadSpecVec = std::vector<std::pair<ReadProps, size_t>>;

ReadVec ReadGen::generate_reads (ReadSpecVec rsv) {

    // NOTE in some kind of weird intermediate state right now
    // with this API
    for (const auto& [spec, n] : rsv) {
        // create partially set bam1_t with templated properties
        // as set by read props
        // (TODO, placeholder below)
        std::string qseq = "AAAA";
        const auto ncigop  = 0;
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

        ReadVec out(n);
        for (auto& b : out) {
          b = bam_init1();  // mandatory
          // bam_set1 (
          //   &b
          //   b_tmpl.l_qname...
          //   next_qname()
          //   and so on
          // )
        }

        return out;
    }

}
