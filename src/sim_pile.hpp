#include <algorithm>
#include <cassert>
#include <cstddef>
#include <format>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <htslib/sam.h>

struct pileup_ev_s {
    size_t A, T, C, G;     // exclusive events
    // size_t del;  // exclusive, not yet implemented
    // size_t fdel, fins, head, tail  // auxillary events
    size_t nreads;

    pileup_ev_s() = delete;

    pileup_ev_s (
        size_t nA,
        size_t nC,
        size_t nG,
        size_t nT
    )
        : A (nA),
          C (nC),
          G (nG),
          T (nT),
          nreads (A + G + C + T) {}
};

struct pileup_props_basic {
    size_t      nread;
    size_t      read_len;
    pileup_ev_s ev;
    std::string ref;

    pileup_props_basic() = delete;

    pileup_props_basic (
        pileup_ev_s ev,
        size_t      read_len,
        std::string ref_template     // NOTE must be double read_len
    )
        : read_len (read_len),
          ev (ev) {
        assert (read_len > 1);     // temporary
        if (ref_template.empty()) {
            ref = std::string (read_len * 2, 'N');
        } else if (ref_template.length() != (read_len * 2)) {
            throw std::runtime_error (
                std::format (
                    "attempt to instantiate pileup_props_basic with "
                    "read length {} not equal to template length {}",
                    read_len,
                    ref_template.length()
                )
            );
        } else {
            ref = ref_template;
        }
    }
};

// NOTE ideally this should be class based for
// lifetime management of the bam1_ts
// NOTE deletions may not be at the end of a query seq
constexpr inline std::vector<bam_pileup1_t>
simulate_pileup (const pileup_props_basic &pr) {
    std::vector<bam_pileup1_t>            out;
    const auto                            ref_pos = pr.read_len - 1;
    std::mt19937                          rng;
    std::uniform_int_distribution<size_t> roll_lmostc (
        0,
        ref_pos
    );     // start coord of reads

    // make shuffle
    std::string evs;     // used as a vector
    evs.append (pr.ev.A, 'A');
    evs.append (pr.ev.T, 'T');
    evs.append (pr.ev.C, 'C');
    evs.append (pr.ev.G, 'G');
    std::shuffle (begin (evs), end (evs), rng);

    while (!evs.empty()) {
        const auto lmc  = roll_lmostc (rng);
        const auto qpos = ref_pos - lmc;
        auto       qseq = pr.ref.substr (lmc, pr.read_len);
        qseq.replace (qpos, 1, 1, evs.back());
        evs.pop_back();

        const auto ncigop  = 1;
        const auto cs      = std::format ("{}M", qseq.length());
        auto       b       = bam_init1();
        uint32_t  *cig_buf = static_cast<uint32_t *> (
            malloc (ncigop * sizeof (uint32_t))
        );
        size_t buf_alloc;
        sam_parse_cigar (cs.c_str(), NULL, &cig_buf, &buf_alloc);
        // bam_parse_cigar (cs.c_str(), NULL, b);
        bam_set1 (
            b,
            0,
            NULL,
            0,
            0,
            lmc,
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

        bam_pileup1_t p1;     // NOTE placeholders abound
        p1.b          = b;
        p1.qpos       = static_cast<int32_t> (qpos);
        p1.indel      = 0;
        p1.level      = 0;
        p1.is_del     = 0;
        p1.is_head    = (qpos == 0);
        p1.is_tail    = (qpos == (qseq.size() - 1));
        p1.is_refskip = 0;
        p1.aux        = 0;
        p1.cigar_ind  = 0;
        // NOTE cd not set
        out.push_back (p1);
    }

    return out;
};
