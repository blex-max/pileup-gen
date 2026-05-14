#include "generate-pileup.hpp"

#include <cstdint>
#include <format>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <string>

#include "plog/Log.h"

#include "read-ops.hpp"
#include "util.hpp"


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
  const auto ref_span = pp.refseq.size();
  return (
    validate (coord)
    && ref_span == span (coord)
    && pp.readlen <= ref_span
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
(const PileupParams& pileup_pars, std::span<const std::pair<size_t, PileupReadSet>> sets)
{
    size_t nsum_reads = 0;
    for (const auto& [nset_reads, _] : sets) {
      nsum_reads += nset_reads;
    }

    PLOGD << std::format ("requested generation of {} total reads", nsum_reads);

    // mem arenas
    PileupData out {
      .b1arr = Bam1Array (new bam1_t[nsum_reads]{}, Bam1ArrayDeleter{nsum_reads}),
      .p1arr = Pileup1Array (new bam_pileup1_t[nsum_reads]{}),
      .nread = nsum_reads
    };

    const auto ref_seq = pileup_pars.refseq;
    const auto pileup_gstart = pileup_pars.coord.gstart;
    const auto pileup_gpos = pileup_pars.coord.gpos;
    const auto pileup_tid = pileup_pars.coord.tid;
    const auto read_len = pileup_pars.readlen;

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

        const auto qpos = set_spec.qpos();
        const auto read_gstart = pileup_gpos - qpos;

        readops::ReadSpec rs{
          .qseq=std::string (genomic_substr (
            pileup_gstart,
            read_gstart,
            read_len,
            ref_seq
          )),  // NOTE: what if the user wants to set specific or empty sequences?
               // well, that wouldn't make sense in this modal context
               // since an empty sequence can't be aligned.
               // NOTE: it doesn't make sense to take the sequence
               // without cigar awareness
          .qqual=set_spec.qual ? (*set_spec.qual)() : std::string (read_len, 37),
          .qname=
            set_spec.qname
            ? (*set_spec.qname)()
            : "read" + std::to_string(i) + "-set" + std::to_string(set_idx),  // NOTE: why would the user ever
                                                                              // care about the qname
                                                                              // in pileup generation mode?
                                                                              // If wanting to generate
                                                                              // corrupted/bad data perhaps,
                                                                              // but it would seem contrived
                                                                              // to be doing so in a pileup
                                                                              // context.
          .qcig={{read_len, readops::cigarcode::match}}, // NOTE: do not be fooled into bothering
                                                         // with overcomplex cigar based generation
                                                         // it will not be rewarded
          .lmost_pos=read_gstart,
          .mate_lmost_pos=read_gstart + 1,  // only relevant for paired end
          .flag=set_spec.flag ? (*set_spec.flag)() : static_cast<uint16_t> (BAM_FREAD1),
          .tid=pileup_tid,
          .mate_tid=set_spec.mate_tid ? (*set_spec.mate_tid)() : pileup_tid,
          .mapq=set_spec.mapq ? (*set_spec.mapq)() : static_cast<uint8_t> (37),
          .aux=set_spec.aux ? (*set_spec.aux)() : std::vector<std::pair<std::string, readops::AuxData>>{}
          // TODO factor out defaults
        };

        // apply perturbation as specified, or no op.
        apply_event (set_spec.event, rs, pileup_gpos);

        readops::set_bam1 (rs, &b1);

        p1 = {
          .b = &b1,
          .qpos = qpos,
          .indel = 0,  // placeholder
          .is_del = 0, // placeholder
          .is_head = (qpos == 0),
          .is_tail = (qpos == read_len - 1),
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
