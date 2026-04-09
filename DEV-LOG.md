# DEV LOG

# blex-max 30/03/26

Started dev log. Performed a significant prune of experimental noodling code
back to a more tractable shape on discussion with Hitham2496. Reflections
copied from a code comment:

"""
For an immediate use case, I want to generate
a pileup with 3 sets of reads, N, M, and R.
Set N reads have variant X and tightly clustered
query positions. Set M and set R reads have equivalently
and uniformly distributed query positions.
Set M reads additionally have variant Y.
Set R reads are identical to reference.

My intuitive sense of the appropriate pattern
is to take a count and a struct of per-property
callables for each set. In other words, A bag of
independently swappable rules for making members of
a set. Note however that with increasing complexity
of the modelled reads some if not all properties
are interdependent upon one another at generation
time. Nevertheless I think the sensible approach
is to make this simple example and then grow from
there rather than trying to overdesign from the outset.

In any case, I'm reasonable sure of the general concept
and design re generation of sets of pileup reads.
A good final design should ensure within reason that
1) sets are modular, 2) it is easy to assemble
overlapping sets, 3) later reuse for common cases is
easy, and 4) users of the library can easily
provide their own generators for use in testing.

The latter point is particularly important. For the
library code, it shifts the question from "what
set of operations is useful to the user" to
"how can we most seamlessly allow user generation
of appropriate data" - which in this specific
context I think is more tractable for this
codebase. And of course it doesn't preclude
providing common shorthands. I think that is
a better approach for this entire project,
not just pileups. i.e. rather than setting up
config machinery for various scenarios like
empty reads, instead make it easy for the
user to python script it themselves.

On that final note, I think centring a CLI
for this kind of small pileup generation was
a misstep. This is library code that should be
used in a scripting language. As such this
subcommand is to be used something like a
script in which we can call the core generation,
without having to commit to figuring out
proper interop at this early stage. On reflection,
I wonder if a CLI may never be the right shape
really for this kind of simulation. It seems
that simulation specification may well be better
off in a high level scripting language, perfomed
by composing functions and helpers.
"""

I'm optimistic that this is a more refined
and product-focused prototype direction.


# blex-max 31/03/26

Did some further stripping, saving an old initial comment
on design philosophy here as it could be useful in the future:

"""
USE CASES:
producing a set of reads aligned to a segment of reference
with "random" noise.
producing a set of reads constituting a pileup at position x
with a known count of each property of interest.
producing a set of reads incorporating target feature/s
e.g. 100 random reads with the motif "ATTTA". Think crispr
quantification.
And more simply, being able to generate a to-spec SAM file
with e.g. a single read pair with feature x, is a common
testing pattern.

FUTURE:
By implementing this functionality well, it should then become
plausible to string them together to generate e.g.
a BAM file of variants to spec

NOTE: it must remain << easier to use this program
than to otherwise assemble the desired data.
Otherwise noone will use it. Beware the endlessly
flexible but incomprehensibly complex config yaml!
"""

The last note is definitely still valid in particular.

Also simplified cigar array generation having found
`bam_cigar_gen ()`.
