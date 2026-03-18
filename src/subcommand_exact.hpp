#include "argparse/argparse.hpp"

#include "generate-reads.hpp"
#include "read-ops.hpp"

inline void setup_exact_parser (argparse::ArgumentParser& args) {
  args.add_argument ("template");
  args.add_argument ("-n")
    .help ("number of reads to produce")
    .default_value (1)
    .nargs(1)
    .scan<'i', int>();  // argparse doesn't support unsigned inputs
}

inline ReadV run_exact (const argparse::ArgumentParser& args) {
  auto template_seq = args.get<std::string> ("template");
  auto n = static_cast<size_t> (args.get<int> ("-n"));

  const auto spec = readops::ReadData{template_seq};

  return generate_reads({{n, spec}});
}
