#pragma once

#include "argparse/argparse.hpp"

#include "generate-reads.hpp"


inline void setup_seq_parser (argparse::ArgumentParser& args) {
  args.add_argument ("ref")
    .help ("string to use as reference");
  args.add_argument ("--read-len")
    .help ("read length")
    .default_value (150)
    .nargs(1)
    .scan<'i', int>();
}

inline ReadV run_seq (const argparse::ArgumentParser& args) {
  auto ref_seq = args.get<std::string> ("ref");
  auto read_len = static_cast<size_t> (args.get<int> ("--read-len"));

  if (read_len > ref_seq.size()) {
    throw std::runtime_error("too long!");
  }

  return {};
}
