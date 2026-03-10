#include <htslib/hts.h>
#include <htslib/sam.h>
#include <iostream>

#include "generate-reads.hpp"
#include "sim_pile.hpp"
#include "util.hpp"
#include "read-ops.hpp"

// USE CASES:
// producing a set of reads constituting a pileup at position x
// with a known count of each property of interest.
//
// producing a set of reads incorporating target feature/s
// e.g. 100 random reads with the motif "ATTTA". Think crispr
// quantification.
// And more simply, being able to generate a to-spec SAM file
// with e.g. a single read pair with feature x, is a common
// testing pattern.
//
// Originally, I envisaged an API for sequential unit testing
// applying edits to a template read, but on reflection
// it is probably conceptually simpler to generate a new one
// based on the same set of arguments modifying equivalently
// to the desired edit operation. In other words
// functionality for creating a single read from template
// arguments.
//
// NOTE: the important point to elucidate for API design
// is how the operations needed for these use cases overlap

// FUTURE:
// By implementing this functionality well, it should then become
// plausible to string them together to generate e.g.
// a BAM file of variants to spec


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


int main(int, char **) {

  // args.add_options(); // todo

  // pileup_ev_s ev(5, 5, 5, 5);
  // pileup_props_basic props(ev, 10, {});

  // auto pile = simulate_pileup(props);

  // // better would be emission in the SAM format
  // // this is temp
  // TODO add this as a subcommand to apb (pileup browser) as `apb print`
  // for a very simple display to terminal
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
  if (const auto rc = sam_hdr_write(hfp, hdr);
      rc < 0) {
    std::cerr << "error writing hdr";
    return 1;
  };

  size_t n = 10;
  const auto spec = readops::ReadArgs{};

  const auto reads = generate_reads({{n, spec}});
  for (const auto& r : reads) {
    if (sam_write1(hfp, hdr, r) < 0) {
      return 1;
    }
  }

  hts_flush(hfp);
  sam_hdr_destroy(hdr);

  return 0;
}
