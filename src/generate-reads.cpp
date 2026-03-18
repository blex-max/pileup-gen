#include <cassert>
#include <htslib/sam.h>
#include <random>
#include <map>

#include "generate-reads.hpp"
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


/*
NOTE: this function is presently retained for
insight into the development of the shape of the program
so the commentary can be written up later
*/
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
    - Chose a sample for that read to have been sequenced from
    - Based on the reference consumption of the edit script,
      choose a starting ref/source position
    - For each simulated read, create a segment edit script
      (basically, an expanded cigar)
      (the edits here will mostly be simulating sequencer error
      since the samples created in the prior step will impart mutation events)
    - Materialise the read

    Pros:
    - Clean separation of sample-level variation and read-level artefacts.
        ^ this is the big one which was tricky while noodling
    - Mutations handled naturally via concrete sample sequence.
    - Segment edit scripts provide a uniform consume/emit model.
    - Avoids coordinate confusion.
    - Extensible to mixtures, haplotypes, and richer error models.

    Cons:
    - Multiple layers add complexity.
    - Risk of over-engineering for current scope.
    - Requires right flank (see below)


    In terms of implementation, the obvious next step is that edit
    script generation can be separated from mutation event generation,
    and creation of fuzzy reads with sequencing error (or any error model,
    I'd suggest a callback here) can be created with a fixed template sequence.

    Note that edit script generation will be context dependent.
    So you need a right flank and/or can only start from ref.size() - read_len.
    There's no get out of jail, that's just how it is!

    Ignore soft clipping for now - it is an aligner feature if sequences
    are highly divergent

    Sequence error events (for model I suppose):
    copy run
    substitution run
    insertion run
    deletion run
    optional untemplated prefix/suffix run (primer/barcode contamination)
    */

}

// presently ignoring the "edit script"
// intermediate representation of
// a read and immediately materialising.
// Can return to that concept later on
// I think.
// I also expect that there will be
// a callback involved for read gen
enum class ReadRegime : size_t {
    none,
    match,
    mismatch
};
struct FromTemplateParams {
    size_t read_len;
};
bam1_t* from_template_sequence (
    std::string_view template_seq,
    std::mt19937& rng,
    const FromTemplateParams& tp
) {
    assert (template_seq.size() >= tp.read_len);
    std::uniform_real_distribution<double> prob_gen (0.0, 1.0);
    std::uniform_int_distribution<size_t> rpos_gen (0, template_seq.size() - tp.read_len);
    const auto rpos_start = rpos_gen (rng);

    // then apply context dependent model
    // FSM
    using RR = ReadRegime;
    // in the future this can be fancier.
    // There are two non-exclusive directions
    // in which to think about this.
    // A regime "op" can be a true operation
    // script describing what happens as you
    // progress through the read - match, mismatch,
    // insertion, deletion.
    // Or a it can be a higher level abstraction
    // describing the state that the model
    // should currently be in based on context
    // e.g. "we are in a homopolymer run"
    // and then it can emit one or more
    // types of specific edit operations.
    // If this machinery ultimately takes
    // callbacks/objects, then both (or more)
    // approaches can be used as long as they
    // return the same result
    using RegimeOp = std::pair<RR, size_t>;
    using RegimeScript = std::vector<RegimeOp>;
    RegimeScript rs;
    RegimeOp curr_op{RR::match, 0};
    const auto finalise_op = [&rs, &curr_op] () {
        rs.emplace_back(curr_op.first, curr_op.second);
        curr_op.second = 0;  // reset len
    };
    for (size_t qpos_counter = 0; qpos_counter < tp.read_len; ++qpos_counter) {
        ++curr_op.second;
        switch (curr_op.first) {
            case (RR::match):
                if (prob_gen (rng) < 0.2) {
                    finalise_op();
                    curr_op.first = RR::mismatch;
                }
                break;
            case (RR::mismatch):
                if (prob_gen (rng) < 0.9) {
                    finalise_op();
                    curr_op.first = RR::match;
                }
                break;
            default:
                throw std::runtime_error ("unknown regime");
        }
    }
    finalise_op();  // finalise last

    // materialise read
    // very simple for now
    readops::ReadData rd;
    auto rpos_cursor = rpos_start;
    for (const auto rop : rs) {
        switch (rop.first) {
            case (RR::match):
                rd.qseq.append(template_seq.substr(rpos_cursor, rop.second));
                rpos_cursor += rop.second;
                break;
            case (RR::mismatch):
                rd.qseq.append(rop.second, 'X');
                rpos_cursor += rop.second;
                break;
            default:
                throw std::runtime_error ("unknown regime");
        }
    }

    return readops::create_bam1(rd);
}

