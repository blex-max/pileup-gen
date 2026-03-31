#include "generate-pileup.hpp"
#include "read-ops.hpp"
#include "util.hpp"
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <string>


size_t span
(const PileupCoordinates& pc) {
  return static_cast<size_t> (pc.gend - pc.gstart);
}
bool validate
(const PileupCoordinates& pc)
{
  // half open 0-based coordinates
  return (
    pc.gstart <= 0
    && pc.gstart < pc.gend
    && pc.gpos >= pc.gstart
    && pc.gpos < pc.gend
    && pc.tid >= 0
  );
}

bool validate
(const PileupParams& pp)
{
  const auto& coord = pp.coord;
  const auto ref_span = pp.ref_region.size();
  return (
    validate (coord)
    && ref_span == span (coord)
    && pp.read_len <= ref_span
  );
}

// TODO/CRITICAL not cigar aware!
// TODO/CRITICAL snp only
void apply_event
(const EventSpec& event, readops::ReadSpec& read, hts_pos_t event_gpos)
{
  if (event.base == BaseEvents::ref) {
    return;
  }

  auto& qseq = read.qseq;
  const auto read_len = qseq.size();
  assert (event_gpos >= read.lmost_pos);
  const auto qpos =
    static_cast <size_t> (event_gpos - read.lmost_pos);
  assert (static_cast<size_t> (qpos) < read_len);

  if (event.base == BaseEvents::del) {
    return;
  } else {
    qseq[qpos] = seq_nt16_str[event.base];
  }
};


// TODO flesh out cigar awareness, application of flags, and anything else
// requiste for the result to properly represent the desired pileup
// explicit specification
PileupData generate_pileup
(const PileupParams& pileup_pars, std::span<const std::pair<size_t, PileupReadSet>> sets, std::mt19937& rng)
{
    size_t nsum_reads = 0;
    for (const auto& [nset_reads, _] : sets) {
      nsum_reads += nset_reads;
    }

    // mem arenas
    PileupData out {
      .b1arr = std::unique_ptr<bam1_t[]> (new bam1_t[nsum_reads]{}),
      .p1arr = std::unique_ptr<bam_pileup1_t[]> (new bam_pileup1_t[nsum_reads]{}),
      .nread = nsum_reads
    };

    const auto ref_seq = pileup_pars.ref_region;
    const auto pileup_gstart = pileup_pars.coord.gstart;
    const auto pileup_gpos = pileup_pars.coord.gpos;
    const auto pileup_tid = pileup_pars.coord.tid;
    const auto read_len = pileup_pars.read_len;

    size_t mem_block_i = 0;
    size_t set_idx = 0;  // name param would be better
    for (const auto& [nread_ev, set_spec] : sets) {
      const auto mem_block_end = nread_ev + mem_block_i;

      for (size_t i=mem_block_i; i < mem_block_end; ++i) {
        auto& b1 = out.b1arr[i];
        auto& p1 = out.p1arr[i];

        /* materialise read according to set spec */

        /*
          NOTE: the create-then-apply flow here is a compromise
          for simplicity and speed. It may be wholly sufficient,
          or it may need revisiting.
          NOTE: can apply sequencing model at creation
        */

        const auto qpos = set_spec.qpos_cb(rng);
        const auto read_gstart = pileup_gpos - qpos;

        // materialised string not string_view
        // since we may then edit it in place.
        readops::ReadSpec rs{
          .qseq=std::string (genomic_substr (
            pileup_gstart,
            read_gstart,
            read_len,
            ref_seq
          )),
          .qqual=std::string (read_len, 37),
          .qname="read" + std::to_string(i) + "-set" + std::to_string(set_idx),  // good candiate/example for callbacks needing
                                              // a context obj of generation up to this point
                                              // + params
          .qcig={{read_len, readops::cigarcode::match}},
          .lmost_pos=read_gstart,
          .mate_lmost_pos=-1,
          .flag=BAM_FREAD1,
          .tid=pileup_tid,
          .mate_tid=-1,
          .mapq=37  // well mapped
        };

        // apply perturbation as specified, or no op.
        apply_event (set_spec.event, rs, pileup_gpos);

        readops::set_bam1 (rs, &b1);

        p1 = {
          .b = &b1,
          .qpos = qpos,
          .indel = 0,
          .is_del = 0,
          .is_head = (qpos == 0),
          .is_tail = (qpos == read_len - 1),
          .is_refskip = 0,
          // NOTE: incomplete initialisation.
          // Likely to factor out to a function later
          // to initialise based on bam1_t + args.
        };
      }

      mem_block_i = mem_block_end;  // move arena ptr
      ++set_idx;

    }

    return out;
};
