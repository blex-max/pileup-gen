#include <cassert>
#include <htslib/sam.h>
#include <random>

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


// trying to generate just one fuzzy read right now
// and build up to multiple.
// I'm imagining perhaps a while loop for each read of diminishing probability
// of adding more events from a randomly generated set of events
// against the reference.
// Also n.b. the ref -> template -> read event propagation idea
ReadV gen_fuzzy_read(std::string_view ref, size_t read_len, std::mt19937& rng) {
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


}

