#pragma once

#include <cassert>
#include <cstddef>
#include <functional>
#include <memory>
#include <span>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "read-ops.hpp"

/* see main.cpp for commentary on
the higher level design and use of this functionality */


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
  std::string ins{};  // .size() == .indel
};


struct PileupCoordinates {
  hts_pos_t gstart;
  hts_pos_t gend;
  hts_pos_t gpos;
  int32_t tid;
};
size_t span (const PileupCoordinates& pc);
bool validate (const PileupCoordinates& pc);


// NOTE WIP
struct PileupParams {
  PileupCoordinates coord;
  std::string_view refseq;
  uint16_t readlen;            // must be <= (ref_region.size() / 2) - 1
};
bool validate (const PileupParams& pp);


// Manifest describing a set
// of reads found in the total
// pileup. Used to materialise
// reads of that set.
struct PileupReadSet {
  EventSpec event;
  // how would these introspect the existing context when used via bindings
  std::function<uint16_t()> qpos_cb;  // callback generating a
                                      // query position from a distribution
                                      // (or otherwise).
  // further properties TODO
  // std::function<std::map<readops::AuxTag, readops::AuxData>()> tag_cb;
};


// mutate a read manifest
// to apply a specified pileup event
void apply_event
(const EventSpec& event, readops::ReadSpec& read, hts_pos_t event_gpos);


// generation output
// functor for freeing mem
struct Bam1ArrayDeleter {
  size_t n;
  void operator()(bam1_t* arr) const {
    for (size_t i = 0; i < n; ++i) {
      auto b = arr[i];
      b.mempolicy = BAM_USER_OWNS_STRUCT;  // don't free structs
      bam_destroy1 (&b);
    }
    delete[] arr;  // free structs
  }
};
using Bam1Array = std::unique_ptr<bam1_t[], Bam1ArrayDeleter>;
using Pileup1Array = std::unique_ptr<bam_pileup1_t[]>;  // does not require custom deleter
struct PileupData {
  // NOTE: destruction will be in reverse order.
  // NOTE: make_unique default-initalises, prefer new T[n]{}
  Bam1Array b1arr;
  Pileup1Array p1arr;
  size_t nread;  // must be set
};
PileupData generate_pileup
(const PileupParams& pileup_pars,
 std::span<const std::pair<size_t, PileupReadSet>> sets,
 PileupReadSet& shared);


