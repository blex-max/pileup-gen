#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <span>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "read-ops.hpp"


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
  // state-free callbacks w.r.t.
  // cpp engine. These may be stateful
  // closures on the python side.
  // Specify an exact set per read profile required,
  // rather than requiring any complex context introspection
  // and stochastics.
  std::function<uint16_t()> qpos;  // callback generating a
                                   // query position from a distribution
                                   // (or otherwise).
  // further properties TODO
  std::optional<std::function<uint16_t()>> flag;
  std::optional<std::function<uint8_t()>> mapq;
  std::optional<std::function<std::string()>> qname; // read{i}-set{j}-{set_name} if not provided
  std::optional<std::function<int32_t()>> mate_tid;  // if not provided, same as pileup tid
                                                     // NOTE: perhaps for user ease this should be
                                                     // a string return which we convert internally
                                                     // to a valid TID
  std::optional<std::function<std::string()>> qual;  // NOTE: AH! this will need to know read len...
                                                     // but again maybe that can be on the python side...
  std::optional<std::function<std::vector<std::pair<std::string, readops::AuxData>>()>> aux;
  std::string set_name;  // optional
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
  Bam1Array b1arr;  // storage backing for pileup1 array
  Pileup1Array p1arr;  // Note this could also be a vector of pointers...
  size_t nread;  // must be set
};
PileupData generate_pileup
(const PileupParams& pileup_pars, std::span<const std::pair<size_t, PileupReadSet>> sets);


