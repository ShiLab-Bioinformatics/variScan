# Positional Index Alignment Algorithm

This program aligns short read sequences against a library of thousands of
closely related reference sequences, such as mutated variants of the same
amplicon, barcode, or genomic locus. It does not perform full dynamic
programming alignment with insertions and deletions. Instead, it finds the best
ungapped placement of each read on each reference. For every reference
sequence, the best mapping location is defined as the placement with the
highest number of matched bases. At that location, the program reports the
start position, number of matched bases, and number of mismatches. This model
is well suited to libraries where variants differ mostly by substitutions.

The key observation is that the reference sequences are highly similar and
mostly the same length. A naive implementation would compare every read against
every reference at every possible start position. For `R` references, read
length `L`, and `S` candidate starts, that costs roughly `R * L * S` base
comparisons per read. With thousands of near-identical references, much of that
work is redundant. The algorithm therefore changes the unit of work: instead of
asking
"how well does this read match reference 1, then reference 2, then reference
3?", it asks "which references have the expected base at this aligned
position?" This allows one lookup to update the score of many references.

To make this possible, the program builds a positional inverted index over the
equal-length references. For each reference position and each base `A/C/G/T`,
the index stores two posting lists: references that match the base at that
position, and references that mismatch it. Conceptually, the structure is:

```text
[match or mismatch] -> [base] -> [reference position] -> [reference IDs]
```

During alignment, a read is placed at a candidate offset. Only the overlapping
part of the read and reference is evaluated, which also allows partial
overlaps near the ends of the reference. For each overlapping query base, the
algorithm looks up the corresponding reference position and updates match or
mismatch counts for many references at once. If a queried base is common at
that position, the matching list may contain most of the library, so it is
cheaper to update only the shorter mismatch list and treat all other references
as implicit matches. If the queried base is rare, the program updates the
shorter match list directly. This adaptive choice matters because most
positions are conserved, while a small number of variant positions distinguish
the references.

For example, suppose the library contains 1,000 reference sequences and the
read has base `A` aligned to reference position 50. If 970 references also have
`A` at position 50, updating the matching list would require 970 counter
updates. It is faster to update only the 30 references that do not have `A` and
then count the other 970 references as implicit matches. At another position,
the same read may carry a rare variant base: only 12 references have `T`, while
988 have another base. In that case, the program updates the 12-reference match
list directly. Both choices give the same match/mismatch counts, but each uses
the shorter posting list.

For each candidate start, the program accumulates two scores for every active
reference: the number of aligned bases that match the read and the number that
mismatch it. After all useful bases for that start have been considered, the
program compares the candidate score with the current best score for that same
reference. If the candidate has more matched bases, it becomes the new best
mapping location for that reference. If the matched-base count is tied, the
program uses the smaller mismatch count as the tie-breaker. In this way, each
reference sequence retains its own best start position, matched-base count, and
mismatched-base count for the read.

The second efficiency feature is the early drop-out rule. Before scanning all
candidate offsets in detail, the program performs a pilot pass against the
first reference sequence and chooses its best read start. Because the library
contains very similar sequences, a good start for the first sequence is usually
informative for the rest of the library. The program then scores all indexed
references at that start to obtain a per-reference lower bound: each reference
already has a score that any later candidate start must improve upon to matter.
At later starts, each reference is kept alive only while it can still beat that
bound. After each processed base, the program computes an optimistic upper
bound:

```text
matches already accumulated + remaining unprocessed bases
```

This is optimistic because it assumes every remaining base will match. If this
upper bound is less than the reference's pilot score, no possible outcome for
the current start can improve that reference's best result. The reference is
marked as stopped for the current start. Once all valid references have stopped,
the candidate start is abandoned immediately. Starts whose overlap length is
already too short to beat the pilot-derived bound are skipped entirely. This
rule is safe because it only drops candidates that cannot mathematically catch
up, even under the most favorable remaining sequence.

References with shorter-than-maximum length are handled separately by direct
scanning, because they cannot share the same fixed-position index as the main
equal-length group. For these sequences, the program tries each possible
ungapped start, counts matches and mismatches directly, and keeps the location
with the highest matched-base count. Even in this direct-scanning path, the
comparison exits early once the mismatch count proves that the current
placement cannot exceed the best match count already seen.

The program is also parallelized across reads with worker goroutines. The index
is built once, shared read-only by all workers, and each read is aligned
independently. This works well because reads do not depend on one another after
the reference index has been constructed. Results are written in a compact
binary format using fixed-width integers, reducing output overhead.

Overall, the algorithm is fast because it exploits the structure of the
problem: references are numerous, similar, and positionally comparable. The
inverted index replaces repeated per-reference base comparisons with batched
posting-list updates, while the early drop-out rule prevents hopeless
reference/start combinations from being scored to completion. The result is an
alignment strategy that is especially efficient for dense panels of related
reference sequences, where conventional all-against-all scanning would spend
most of its time rediscovering the same conserved bases.
