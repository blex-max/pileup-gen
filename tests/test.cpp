#include <catch2/catch_test_macros.hpp>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <iostream>

#include "sim_pile.hpp"


SCENARIO ("create edited reads from template read") {
    GIVEN ("A template read bam1_t tr;") {
        auto tr = bam_init1();
        // bam_set1(
        //     tr
        // )
        WHEN ("clone_read(tr, ...) is called with no edits") {
            // const auto cb = from_template(tr);
            THEN ("The read should be identical to the template") {
                // REQUIRE(bam1_cmp(tr, cb))
            }
        }

        WHEN ("clone_read(tr, ...) is called with a simple substitution edit") {
            THEN ("The subsitution is present on the query sequence") {

            }
            THEN ("The read exactly matches the template, other than the substitution") {

            }
        }

        WHEN ("clone_read(tr, ...) is called with a deletion edit") {
            THEN ("The deletion is represented on the query sequence") {

            }
            THEN ("The deletion is represented on the quality sequence") {

            }
            THEN ("The deletion is represented in the cigar operations") {

            }
            THEN ("The read exactly matches the template, other than the substitution") {

            }
        }
    }
}


SCENARIO ("simple simulated pileup") {
    size_t Abases = 6;
    size_t Cbases = 1;
    size_t Gbases = 3;
    size_t Tbases = 0;
    GIVEN (
        std::format (
            "A simple set of exclusive alignment events:\n\t"
            "A:{}, C:{}, G:{}, T:{}\n",
            Abases,
            Cbases,
            Gbases,
            Tbases
        )
    ) {
        pileup_props_basic props ({6, 1, 3, 0}, 10, "");
        WHEN ("simulate_pileup() is called") {
            // generate data
            const std::vector<bam_pileup1_t> reads = simulate_pileup (
                props
            );
            THEN ("A set of reads must be returned matching that set of exclusive events") {
                std::array<size_t, 4> observed_counts{0, 0, 0, 0};
                std::array<size_t, 4> expected_counts{
                    props.ev.A,
                    props.ev.C,
                    props.ev.G,
                    props.ev.T
                };
                // check (and print pileup for display)
                std::cout << "pileup visual:" << std::endl;
                std::cout << props.ref << " (ref)" << std::endl;
                for (const auto &b : reads) {
                    const auto  bseq = bam_get_seq (b.b);
                    std::string qs (b.b->core.pos, ' ');     // offset
                    for (size_t i = 0; i < b.b->core.l_qseq; ++i) {
                        qs.push_back (seq_nt16_str[bam_seqi (bseq, i)]);
                    }
                    observed_counts
                        [seq_nt16_int[bam_seqi (bam_get_seq (b.b), b.qpos)]]++;
                    std::cout << qs << std::endl;
                }
                UNSCOPED_INFO (
                    "testing the number of reads returned "
                    "matches the total number of events requested."
                );
                REQUIRE (reads.size() == props.ev.nreads);
                UNSCOPED_INFO ("testing that each event occurs the number of times requested.");
                REQUIRE (observed_counts == expected_counts);
            }
        }
    }
}


// TODO
// generate new read from template read (the copy/edit approach I use in python) - see github discussion
