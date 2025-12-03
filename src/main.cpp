#include <cxxopts.hpp>
#include <htslib/sam.h>

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
// -> could do this in pysam (might be trickier)

// TODO
// plan out each unit/function/module needed

int main(int argc, char *argv[]) { return 0; }
