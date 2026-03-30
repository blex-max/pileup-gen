#include <plog/Log.h>
#include <plog/Initializers/ConsoleInitializer.h>
#include <plog/Formatters/TxtFormatter.h>

#include <argparse/argparse.hpp>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "sim_pile.hpp"
#include "util.hpp"
#include "subcommand_exact.hpp"
#include "subcommand_seq.hpp"
#include "subcommand_pileup.hpp"

// USE CASES:
// producing a set of reads aligned to a segment of reference
// with "random" noise.
//
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
//
// NOTE: it must remain << easier to use this program
// than to otherwise assemble the desired data.
// Otherwise noone will use it. Beware the endlessly
// flexible but incomprehensibly complex config yaml!


// NOTE: cli shape is demonstrative, not obligatory,
// and will almost certainly need modification
// in response to pratical challenges met
// during implementation
// (as with all code written so far I suppose).
int main (int argc, char** argv) {
  argparse::ArgumentParser cli ("htsgen", "0.0.0");

  argparse::ArgumentParser sub_exact ("exact");
  sub_exact.add_description ("generate reads according to exact spec");
  setup_exact_parser (sub_exact);
  cli.add_subparser (sub_exact);

  argparse::ArgumentParser sub_pileup ("pileup");
  sub_pileup.add_description ("generate reads by specifying a pileup");
  setup_pileup_parser(sub_pileup);
  cli.add_subparser(sub_pileup);
  
  argparse::ArgumentParser sub_seq ("seq");
  sub_seq.add_description ("generate reads via simulation of sequencing");
  setup_seq_parser (sub_seq);
  cli.add_subparser(sub_seq);

  try {
    cli.parse_args(argc, argv);
  }
  catch (const std::exception& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << cli;
    return 1;
  }

  plog::init<plog::TxtFormatter>(plog::debug, plog::streamStdErr);

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

  bam1_t* read_arr;
  size_t nread = 0;
  if (cli.is_subcommand_used(sub_exact)) {
    PLOGD << "exact subcommand called";
    try {
      // reads = run_exact (sub_exact);
    }
    catch (const std::exception& ex) {
      PLOGF << ex.what();
      return 1;
    }
  }
  if (cli.is_subcommand_used(sub_seq)) {
    PLOGD << "seq subcommand called";
    try {
      // reads = run_seq (sub_seq);
    }
    catch (const std::exception& ex) {
      PLOGF << ex.what();
      return 1;
    }
  }
  if (cli.is_subcommand_used(sub_pileup)) {
    PLOGD << "pileup subcommand called";
    try {
      // pileup is statically allocated because
      // the return type PileupData holds
      // unique ptrs, which need to live till end
      // of program. So much for streamlined memory management.
      // Might as well just pass in a preallocated block of
      // max_depth or similar. Don't think about this too hard
      // until the library shape and higher level bindings are
      // available though.
      static auto pileup = run_pileup (sub_pileup);
      read_arr = pileup.b1arr.get();
      nread = pileup.nread;
    }
    catch (const std::exception& ex) {
      PLOGF << ex.what();
      return 1;
    }
  }

  PLOGD << "writing reads";
  for (size_t i = 0; i < nread; ++i) {
    if (sam_write1(hfp, hdr, read_arr + i) < 0) {
      PLOGF << "failed to write";
      return 1;
    }
  }

  hts_flush(hfp);
  hts_close(hfp);
  sam_hdr_destroy(hdr);

  return 0;
}
