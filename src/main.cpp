#include <cstdlib>
#include <fstream>
#include <random>
#include <filesystem>

#include <plog/Log.h>
#include <plog/Initializers/ConsoleInitializer.h>
#include <plog/Formatters/TxtFormatter.h>

#include <argparse/argparse.hpp>

#include <htslib/hts.h>
#include <htslib/sam.h>

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

  On that final note, I think centring a CLI
  for this kind of small pileup generation was
  a misstep. This is library code that should be
  used in a scripting language. As such this
  subcommand is to be used something like a
  script in which we can call the core generation,
  without having to commit to figuring out
  proper interop at this early stage. On reflection,
  I wonder if a CLI may never be the right shape
  really for this kind of simulation. It seems
  that simulation specification may well be better
  off in a high level scripting language, perfomed
  by composing functions and helpers.
*/

namespace fs = std::filesystem;


int main (int argc, char** argv)
{
  plog::init<plog::TxtFormatter> (plog::debug, plog::streamStdErr);
  argparse::ArgumentParser cli ("htsgen", "0.0.0");

  cli.add_argument ("outdir")
    .help ("writeable directory in which to output results")
    .nargs (1);
  cli.add_argument ("--prefix")
    .help ("filename prefix to add to all outputs")
    .nargs (1);

  try {
    cli.parse_args(argc, argv);
  }
  catch (const std::exception& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << cli;
    return EXIT_FAILURE;
  }

  const auto fprefix = cli.get<std::string> ("--prefix");

  fs::path outdir = cli.get<std::string> ("outdir");
  if (fs::exists (outdir) && !fs::is_directory (outdir)) {
    std::cerr << "outdir exists and is not a directory" << std::endl;
    return EXIT_FAILURE;
  }

  fs::create_directory (outdir);  // no-op if exists
  if ((fs::status (outdir).permissions() & fs::perms::owner_write) == fs::perms::none) {
    std::cerr << "outdir not writeable" << std::endl;
    return EXIT_FAILURE;
  }

  fs::path ref_out = outdir / (fprefix + "pileup-ref.fa");
  std::ofstream ref_ostream (ref_out);
  if (ref_ostream.fail()) {
    std::cerr << "failed to open output fasta for writing at " << ref_out;
    return EXIT_FAILURE;
  }

  fs::path sam_out = outdir / (fprefix + "pileup.sam");
  auto sam_handle = hts_open (sam_out.c_str(), "w");
  auto hdr = sam_hdr_init ();
  // NOTE attempting to write a bam1_t
  // to file without having set SQ lines
  // is a segfault if any RNAME
  // is set in the bam1_t
  const std::string tid ("chr1");
  sam_hdr_add_line (hdr, "SQ",
                   "SN", tid.c_str(),
                   "LN", "248956422",
                   NULL);
  if (const auto rc = sam_hdr_write (sam_handle, hdr);
      rc < 0) {
    PLOGF << "error writing hdr";
    return EXIT_FAILURE;
  };

  /* setup pileup description */
  std::mt19937 rng;
  const uint16_t read_len = 50;
  const std::string ref (read_len + read_len, 'G');
  const size_t nreads_alt = 20;
  const size_t nreads_ref = 60;
  PileupParams ppars {
    .coord={
      .gstart=0,
      .gend=static_cast<hts_pos_t> (ref.size()),
      .gpos= read_len - 1,
      .tid=sam_hdr_name2tid(hdr, tid.c_str())
    },
    .refseq=ref,
    .readlen=read_len
  };
  if (!validate (ppars)) {
    PLOGF << "invalid pileup specification";
    return EXIT_FAILURE;
  };

  /* specifications for each set of reads I want to find
     in the pileup
  */
  // TODO output metadata
  std::uniform_int_distribution<uint16_t> broad_ud (0, read_len - 1);
  assert (read_len > 2);
  const uint16_t read_midpoint= (read_len / 2) - 1;
  const auto wobble = static_cast<uint16_t> (ceil (read_len * 0.05));
  std::uniform_int_distribution<uint16_t> clust_ud (
    read_midpoint - wobble , read_midpoint + wobble
  );
  PileupReadSet set_a {
    .event={BaseEvents::A},
    .qpos_cb=[&broad_ud, &rng] () { return broad_ud(rng); }
  };
  // PileupReadSet set_b {
  //   .event={BaseEvents::T},
  //   .qpos_cb=[&clust_ud, &rng] () { return clust_ud(rng); }
  // };
  PileupReadSet set_ref {
    .event={BaseEvents::ref},
    .qpos_cb=[&broad_ud, &rng] () { return broad_ud(rng); }
  };

  /* generate */
  const std::vector<std::pair<size_t, PileupReadSet>> evs
  {
    {nreads_alt, set_a},
    // {nreads_alt, set_b},
    {nreads_ref, set_ref}
  };
  const auto pileup = generate_pileup (ppars, evs);
  const auto read_arr = pileup.b1arr.get();

  PLOGD << "writing reads";
  for (size_t i = 0; i < pileup.nread; ++i) {
    const auto rc = sam_write1(sam_handle, hdr, read_arr + i);
    if (rc < 0) {
      PLOGF << std::format ("failed to write bam1_t, error code {}", rc);
      return EXIT_FAILURE;
    }
  }

  if (ref_ostream.is_open()) {
    PLOGD << "writing reference";
    ref_ostream << ">" << tid << "\n";
    ref_ostream << ppars.refseq << "\n";
    ref_ostream.flush();
    ref_ostream.close();
  }

  hts_flush(sam_handle);
  hts_close(sam_handle);
  sam_hdr_destroy(hdr);

  return EXIT_SUCCESS;
}
