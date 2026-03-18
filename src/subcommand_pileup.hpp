#pragma once

#include "argparse/argparse.hpp"
#include "generate-reads.hpp"

inline void setup_pileup_parser (argparse::ArgumentParser& args) {
  // TODO
}

inline ReadV run_pileup (const argparse::ArgumentParser& args) {
  // TODO
  // pileup_ev_s ev(5, 5, 5, 5);
  // pileup_props_basic props(ev, 10, {});
  // auto pile = simulate_pileup(props);
  // TODO add this as a subcommand to apb (pileup browser) as `apb print`
  // for a very simple display to terminal
  // std::cout << std::format ("{}", props.ref) << "\n";
  // for (const auto r : pile) {
  //   const std::string pad(static_cast<size_t> (r.b->core.pos), ' ');
  //   const auto seq = get_seq(r.b);
  //   std::cout << std::format ("{}{}", pad, seq) << "\n";
  // }

  return {};
}
