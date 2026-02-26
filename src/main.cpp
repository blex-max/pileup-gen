#include <cxxopts.hpp>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <iostream>

#include "sim_pile.hpp"
#include "util.hpp"
#include "read-ops.hpp"

// SKETCH:
// minimal:
// from an input array of event counts at position x
// produce a set of reads which fulfills those counts at pileup postion x

// at CLI:
// from an input tsv/csv of events + config
// + optional bam header to copy
// + optional reference fasta + region spec
// and output depending on mode:
// -> produce a fastq of reads
// -> produce a bam file of unaligned reads
// -> produce a bam file of aligned reads
// Quick shorthands for common use cases


// NOTES:
// -> could provide multiple sparse genomic positions of events
// to produce alignment over region
// -> could provide multiple sparse query positions of events
// to produce specified multi event reads
// -> can apply models to read simulation
// -> could do this in pysam (might be trickier)

// TODO
// plan out each unit/function/module needed

int main(int argc, char *argv[]) {

  cxxopts::Options args("hts-syn", "simulate htslib pileup");

  // args.add_options(); // todo

  // pileup_ev_s ev(5, 5, 5, 5);
  // pileup_props_basic props(ev, 10, {});

  // auto pile = simulate_pileup(props);

  // // better would be emission in the SAM format
  // // this is temp
  // std::cout << std::format ("{}", props.ref) << "\n";
  // for (const auto r : pile) {
  //   const std::string pad(static_cast<size_t> (r.b->core.pos), ' ');
  //   const auto seq = get_seq(r.b);
  //   std::cout << std::format ("{}{}", pad, seq) << "\n";
  // }

  auto hfp = hts_open("-", "w");
  auto hdr = sam_hdr_init();
  // NOTE attempting to write a bam1_t
  // to file without having set SQ lines
  // is a hard segfault if any RNAME
  // is set in the bam1_t
  // sam_hdr_add_line(hdr, "SQ",
  //                  "SN", "chr1",
  //                  "LN", "248956422",
  //                  NULL);    
  sam_hdr_write(hfp, hdr);

  // auto b = bam_init1();

  size_t n = 10;
  const auto spec = readops::ReadArgs{};
  for (size_t i = 0; i < n; ++i) {
    // TODO use smart pointers
    const auto b = readops::create_read(spec);
    // sam_write1(hfp, hdr, b);
    if (sam_write1(hfp, hdr, b) < 0) {
      return 1;
    }
    bam_destroy1(b);
  }

  hts_flush(hfp);

  // bam_destroy1(b);
  sam_hdr_destroy(hdr);

  return 0;
}
