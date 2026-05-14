"""
Walkthrough: generate a small synthetic pileup with two populations of reads.

  - 20 "alt" reads carrying an A at a broadly-distributed query position
  - 60 "ref" reads identical to reference, also broadly distributed

"""

import random
import htsgen
from collections.abc import Callable


# --- specify pileup --- #
READ_LEN = 50
REF_SEQ: str = "G" * (READ_LEN * 2)   # 100 bp all-G reference
TID = 0                       # chr1 is the only contig; numeric id = 0
NREADS_ALT = 20
NREADS_REF = 60

coords = htsgen.PileupCoordinates(
    gstart = 0,
    gend = len(REF_SEQ),  # off by one?
    gpos = READ_LEN - 1,   # position of interest
    tid = TID,
)  # TODO: alternative constructor taking only region length and tid, to infer the rest

ppars = htsgen.PileupParams(
    coordinates = coords,
    refseq = REF_SEQ,
    readlen = READ_LEN,
)  # TODO: wrappr that constructs params and coordinates from minimal info

# --- specify read sets --- #
# Broad: uniform query position across the full read length
broad_qpos = lambda: random.randint(0, READ_LEN - 1)

# Clustered
midpoint = (READ_LEN // 2) - 1
wobble   = int(READ_LEN * 0.05)
clust_qpos = lambda: random.randint(midpoint - wobble, midpoint + wobble)

# all our reads will have this same tag data
def make_aux_tags () -> list[tuple[str, str | int | float]]:
    return [
        ("MC", "50M"),
        ("AS", 100)
    ]

# helper
def PileupReadSet_factory (
    event: htsgen.EventSpec,
    fn_qpos: Callable[[], int]
) -> htsgen.PileupReadSet:
    """
    factory to produce read sets
    with specified properties and
    fixed shared properties.
    """
    return htsgen.PileupReadSet (
        event=event,
        qpos=fn_qpos,
        aux=make_aux_tags
    )

# a broadly distributed variant allele
set_a = PileupReadSet_factory (
    htsgen.EventSpec(htsgen.BaseEvents.A),
    broad_qpos,
)

# a clustered variant allele
set_t = PileupReadSet_factory (
    htsgen.EventSpec(htsgen.BaseEvents.T),
    clust_qpos,
)

# broadly distributed reference allele
set_ref = PileupReadSet_factory (
    htsgen.EventSpec(htsgen.BaseEvents.ref),
    broad_qpos,
)

# --- generate --- #
pileup = htsgen.generate_pileup(
    ppars,
    [(NREADS_ALT, set_a),
    (NREADS_ALT, set_t),
    (NREADS_REF, set_ref)]
)

# --- use result --- #
# Do work with the generated pileup.
# PileupData is iterable; each entry is a PileupEntry,
# which is a thin wrapper around bam_pileup1_t
#
# e.g:
# print(f"{'base':<6} {'qpos':<6} {'gstart':<8} {'qual':<6} is_del")
# print("-" * 38)
# for entry in pileup:
#     print(f"{entry.base:<6} {entry.qpos:<6} {entry.gstart:<8} {entry.base_qual:<6} {entry.is_del}")

# and/or write out
print ("nreads: {}".format(pileup.nread));
print ("writing generated data")
with open("pileup-ref.fa", "w") as fh:
    _ = fh.write(">chr1")
    _ = fh.write(REF_SEQ)
_ = htsgen.write_pileup(pileup, "pileup.sam");

