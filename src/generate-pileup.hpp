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

/*
  For an immediate use case, I want to generate
  a pileup with 3 sets of reads, N, M, and R.
  Set N reads have variant X and tightly clustered
  query positions. Set M and set R reads have equivalently
  and uniformly distributed query positions.
  Set M reads additionally have variant Y.
  Set R reads are identical to reference.

  My intuitive sense of the appropriate pattern
  is to take a count and a struct of per-property
  callables for each set. In other words, A bag of
  independently swappable rules for making members of
  a set. Note however that with increasing complexity
  of the modelled reads some if not all properties
  are interdependent upon one another at generation
  time. Nevertheless I think the sensible approach
  is to make this simple example and then grow from
  there rather than trying to overdesign from the outset.

  In any case, I'm reasonable sure of the general concept
  and design re generation of sets of pileup reads.
  A good final design should ensure within reason that
  1) sets are modular, 2) it is easy to assemble
  overlapping sets, 3) later reuse for common cases is
  easy, and 4) users of the library can easily
  provide their own generators for use in testing.

  The latter point is particularly important. For the
  library code, it shifts the question from "what
  set of operations is useful to the user" to
  "how can we most seamlessly allow user generation
  of appropriate data" - which in this specific
  context I think is more tractable for this
  codebase. And of course it doesn't preclude
  providing common shorthands.
*/


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
struct PileupSharedParams {
  PileupCoordinates coord;
  // can provide any single, and the others will be derived
  std::string_view ref_region;  // ref pos is center, therefore N must be odd
  uint16_t read_len;         // <= (ref_region.size() + 1) / 2
};


// output
struct PileupData {
  // NOTE: destruction will be in reverse order.
  // make_unique default-initalises, prefer new T[n]{}
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


// explicit specification
inline PileupData generate_pileup
(const PileupSharedParams& pileup_pars, std::span<const std::pair<size_t, PileupReadSet>> sets, std::mt19937& rng) {
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

    size_t mem_block_i= 0;
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

      mem_block_i = mem_block_end;

    }

    return out;
};


// random
// PileupVec
// generate (double pt_non_ref, size_t nreads);

// // percentage spec
// PileupVec
// generate (const std::vector<std::pair<double, EventSpec>> &events, size_t nreads);

