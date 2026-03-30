#pragma once

#include <cstdint>
#include <random>

#include "argparse/argparse.hpp"
#include "generate-pileup.hpp"


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
  providing common shorthands. I think that is
  a better approach for this entire project,
  not just pileups. i.e. rather than setting up
  config machinery for various scenarios like
  empty reads, instead make it easy for the
  user to python script it themselves.
*/


inline void setup_pileup_parser (argparse::ArgumentParser&) {
  // TODO
}


inline PileupData run_pileup
(const argparse::ArgumentParser&)
{
  /* basic setup */
  std::mt19937 rng;
  const uint16_t read_len = 50;
  const std::string ref (read_len + read_len, 'G');
  const size_t nreads = 10;
  PileupSharedParams ppars {
    .coord={
      .gstart=0,
      .gend=static_cast<hts_pos_t> (ref.size()),
      .gpos= read_len - 1
    },
    .ref_region=ref,
    .read_len=read_len
  };
  // TODO validate_pileup_pars (ppars);

  /* specifications for each set of reads I want to find
     in the pileup
  */
  std::uniform_int_distribution<uint16_t> broad_ud (0, read_len);
  assert (read_len > 2);
  const uint16_t read_midpoint= (read_len / 2) - 1;
  const auto wobble = static_cast<uint16_t> (ceil (read_len * 0.05));
  std::uniform_int_distribution<uint16_t> clust_ud (
    read_midpoint - wobble , read_midpoint + wobble
  );
  PileupReadSet set_a {
    .event={BaseEvents::A},
    .qpos_cb=[&broad_ud] (std::mt19937& rng) { return broad_ud(rng); }
  };
  PileupReadSet set_b {
    .event={BaseEvents::T},
    .qpos_cb=[&clust_ud] (std::mt19937& rng) { return clust_ud(rng); }
  };
  PileupReadSet set_ref {
    .event={BaseEvents::ref},
    .qpos_cb=[&broad_ud] (std::mt19937& rng) { return broad_ud(rng); }
  };

  /* generate */
  const std::vector<std::pair<size_t, PileupReadSet>> evs
  {
    {nreads, set_a},
    {nreads, set_b},
    {nreads, set_ref}
  };
  // NOTE: returned reads are somewhat incomplete, e.g unalgined
  return generate_pileup (ppars, evs, rng);
}
