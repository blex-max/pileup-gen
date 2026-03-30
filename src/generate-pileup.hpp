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

/* see "subcommand-pileup.hpp" for commentary on
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


void apply_event
(const EventSpec& event, readops::ReadSpec& read, hts_pos_t event_gpos);


PileupData generate_pileup
(const PileupParams& pileup_pars, std::span<const std::pair<size_t, PileupReadSet>> sets, std::mt19937& rng);


