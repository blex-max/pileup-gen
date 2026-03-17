#include <cassert>
#include <htslib/sam.h>
#include <random>
#include <map>

#include "generate-reads.hpp"
#include "generate-pileup.hpp"  // must factor out eventspec!!
#include "read-ops.hpp"


ReadV generate_reads (const ReadArgV& rav) {
    ReadV out;
    for (const auto& [n, args] : rav) {
        for (size_t i = 0; i < n; ++i) {
            out.push_back (readops::create_bam1(args));
        }
    }
    return out;
}


enum MutationType: size_t {
    transition,
    ins,
    del
};
struct MutationEvent {
    MutationType t;
    size_t len;
    float intensity=1.0;  // 0-1 likelihood
};

// trying to generate just one fuzzy read right now
// and build up to multiple.
// I'm imagining perhaps a while loop for each read of diminishing probability
// of adding more events from a randomly generated set of events
// against the reference.
// Also n.b. the ref -> template -> read event propagation idea
ReadV gen_fuzzy_reads
(
    size_t n,
    std::string_view ref,
    size_t read_len,
    float event_intensity,  // 0 - 1, likelihood per base (there are many ways one could do this)
    std::mt19937& rng
) {
    assert (ref.size() >= read_len);

    // A read can be described as a series of events
    // this manifest can then be iterated through and the read
    // materialised.
    // An empty optional is "no event" - matching reference, normal looking base.
    // The point of doing it this way is that it allows creating a full
    // read agnostic of starting coordinate. If we materialise immediately
    // off of reference it becomes tricky to consider how soft-clipping
    // and other non reference-consuming events fit in.
    // std::vector<std::optional<EventSpec>> read_manifest(read_len);
    // The above is kind of like a more abstracted version of a cigar string.
    // Being more abstracted, it's harder to think about. For a first pass,
    // perhaps see below:

    // generate a cigar, then generate a read
    // [const auto cig, const auto ref_consumed] = gen_cigar (read_len, rng)
    // const auto lc = uniform_int_distribution<size_t> (0, ref.size() - ref_consumed)
    // const auto qseq = seq_from_cigar(cig, ref.substr(lc, ref_consumed))  // what will this do about insertions?
    // really this should be filling out a ReadData struct
    // const auto qqual = gen_qual(cig, qseq, ref) // reality doesn't cascade so smoothly..., which is the appeal
    // of the higher level event based absraction I suppose
    //
    // Also, most reads will not have independent large differences, which is the appeal of the
    // sample -> template/fragment -> reads approach. Or in other words, it seems
    // sensible to generate a fixed set of divergences from the reference, and then
    // apply a randomly selected set of those, +- sequencing error, to generate reads.

    // gen_events

    const auto n_ev =
        static_cast<size_t> (
            ceil (event_intensity
                  * static_cast<float> (ref.size())));
    std::map<size_t, MutationEvent> sample_ev;  // keyed by genomic position, ascending
    auto position_gen = std::uniform_int_distribution<size_t> (0, ref.size());
    auto dice_gen = std::uniform_int_distribution<size_t> (0, 100);
    for (size_t i = 0; i < n_ev; ++i) {
        const auto epos = position_gen(rng);
        // keep it simple for now
        sample_ev.try_emplace(epos, MutationType::transition, 1);
        // TODO e.g.
        // const auto dice_roll = dice_gen(rng);
        // if (dice_roll > 30) {
        //     sample_ev.try_emplace(epos, MutationType::transition, 1);
        // } else if (dice_roll > 15) {
        //     sample_ev.try_emplace(epos, MutationType::del, 3);
        // } else {
        //     sample_ev.try_emplace(epos, MutationType::ins, 3);
        // }
    }

    // TODO this approach does NOT
    // account for non ref consuming events (see above comments).
    auto qpos_gen = std::uniform_int_distribution<size_t> (0, ref.size() - read_len);
    for (size_t i = 0; i < n; ++i) {
        const auto qpos = qpos_gen(rng);
        // get overlapping events
        auto it_low = sample_ev.lower_bound(qpos);
        auto it_hi = sample_ev.upper_bound(qpos + read_len);
        for (auto it = it_low; it != it_hi; ++it) {
            // apply event
        }
        // readops::create_bam1(const ReadData &ra)
    }

    /*
    TODO/NEXT
    The above was noodling.
    Noodling has led to the following conclusion/
    rough shape for this process:
    - Create a list of sample mutation events
    - Create one or more sample sequences, edited from reference,
      applying some or all of those events
    - For each simulated read, create a segment edit script
      (basically, an expanded cigar)
      (the edits here will mostly be simulating sequencer error
      since the samples created in the prior step will impart mutation events later)
    - Chose a sample for that read to have been sequenced from
    - Based on the reference consumption of the edit script,
      choose a starting ref/source position
    - Materialise the read

    Pros:
    - Clean separation of sample-level variation and read-level artefacts.
    - Indels and clipping handled naturally via concrete sample sequence.
    - Segment edit scripts provide a uniform consume/emit model.
    - Avoids coordinate confusion.
    - Extensible to mixtures, haplotypes, and richer error models.

    Cons:
    - Multiple layers add complexity.
    - Risk of over-engineering for current scope.
    */

}

