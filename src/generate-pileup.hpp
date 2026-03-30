#pragma once

#include <cassert>
#include <cstddef>
#include <functional>
#include <memory>
#include <random>
#include <span>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "read-ops.hpp"
#include "util.hpp"


// maps to ambiguity codes
// plus del.
// probably better to use
// cigar ops and separate
// base/mismatch table
enum BaseEvents : size_t {
  ref,
  A=1,
  C=2,
  G=4,
  T=8,
  N=16,
  del=17
};


// concrete event
struct EventSpec {
  BaseEvents base;
  int indel=0;        // <0 del, >0 ins, length indel
  std::string ins{};  // empty is random, mismatched length is err
};


struct PileupCoordinates {
  hts_pos_t gstart;
  hts_pos_t gend;
  hts_pos_t gpos;
};


// NOTE WIP
struct PileupParams {
  PileupCoordinates coord;
  std::string_view ref_region;  // ref pos is center, therefore N must be odd
  uint16_t read_len;            // must be <= (ref_region.size() / 2) - 1
};


// output
struct PileupData {
  // NOTE: destruction will be in reverse order.
  // NOTE: make_unique default-initalises, prefer new T[n]{}
  std::unique_ptr<bam1_t[]> b1arr;
  std::unique_ptr<bam_pileup1_t[]> p1arr;
  size_t nread;  // must be set
};


struct PileupReadSet {
  EventSpec event;
  std::function<uint16_t(std::mt19937&)> qpos_cb;  // callback
};


// TODO/CRITICAL not cigar aware!
// TODO/CRITICAL snp only
inline void apply_event
(const EventSpec& event, readops::ReadSpec& read, hts_pos_t event_gpos) {
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
inline PileupData generate_pileup
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
    const auto read_len = pileup_pars.read_len;

    size_t mem_block_i = 0;
    for (const auto& [nread_ev, set_spec] : sets) {
      const auto mem_block_end = nread_ev + mem_block_i;

      for (size_t i=mem_block_i; i < mem_block_end; ++i) {
        auto& b1 = out.b1arr[i];
        auto& p1 = out.p1arr[i];

        const auto qpos = set_spec.qpos_cb(rng);
        const auto rstart = pileup_gpos - qpos;

        // materialised string not string_view
        // since we may then edit it in place.
        readops::ReadSpec rs{
          .qseq=std::string (genomic_substr (
            pileup_gstart,
            rstart,
            read_len,
            ref_seq
          )),
          .qqual=std::string (read_len, 37),
          .qname={},
          .qcig={},
          .lmost_pos=rstart,
          .mate_lmost_pos=0,
          .flag=BAM_FUNMAP,
          .tid=0,
          .mate_tid=0,
          .mapq=0
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

    }

    return out;
};


// random
// PileupVec
// generate (double pt_non_ref, size_t nreads);

// // percentage spec
// PileupVec
// generate (const std::vector<std::pair<double, EventSpec>> &events, size_t nreads);

