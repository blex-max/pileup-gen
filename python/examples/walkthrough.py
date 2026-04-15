"""
Walkthrough: generate a small synthetic pileup with two populations of reads.

  - 20 "alt" reads carrying an A at a broadly-distributed query position
  - 60 "ref" reads identical to reference, also broadly distributed

"""

import random
import htsgen


# --- specify pileup --- #
READ_LEN = 50
REF_SEQ: str = "G" * (READ_LEN * 2)   # 100 bp all-G reference
TID = 0                       # chr1 is the only contig; numeric id = 0
NREADS_ALT = 20
NREADS_REF = 60

coords = htsgen.PileupCoordinates(
    gstart = 0,
    gend = len(REF_SEQ),
    gpos = READ_LEN - 1,   # position of interest
    tid = TID,
)

ppars = htsgen.PileupParams(
    coordinates = coords,
    refseq = REF_SEQ,
    readlen = READ_LEN,
)

# --- specify read sets --- #
# Broad: uniform query position across the full read length
broad_qpos = lambda: random.randint(0, READ_LEN - 1)

# Clustered (mirrors commented-out set_b in main.cpp):
# midpoint = (READ_LEN // 2) - 1
# wobble   = int(READ_LEN * 0.05)
# clust_qpos = lambda: random.randint(midpoint - wobble, midpoint + wobble)

# The only supported "interesting"
# functionality right now is the ability to
# specify a callback for query position
set_a = htsgen.PileupReadSet(
    event = htsgen.EventSpec(htsgen.BaseEvents.A),
    qpos_cb = broad_qpos,
)

set_ref = htsgen.PileupReadSet(
    event = htsgen.EventSpec(htsgen.BaseEvents.ref),
    qpos_cb = broad_qpos,
)

# --- generate --- #
pileup = htsgen.generate_pileup(ppars, [
    (NREADS_ALT, set_a),
    (NREADS_REF, set_ref),
])

# --- use result --- #
# Do work with the generated pileup.
# PileupData is iterable; each entry is a PileupEntry,
# which is a thin wrapper around bam_pileup1_t
print(f"{'base':<6} {'qpos':<6} {'gstart':<8} {'qual':<6} is_del")
print("-" * 38)
for entry in pileup:
    print(f"{entry.base:<6} {entry.qpos:<6} {entry.gstart:<8} {entry.base_qual:<6} {entry.is_del}")

