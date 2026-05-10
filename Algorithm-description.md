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

## Data Structure
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
[match or mismatch] -> [reference position] -> [base] -> [reference IDs]
```

## Calculation of Match/Mismatch Counts against All Reference Sequences Simutaniously
During alignment, a read is placed at a candidate offset. Only the overlapping
part of the read and reference is evaluated, which also allows partial
overlaps near the ends of the reference. For each overlapping query base, the
algorithm looks up the corresponding reference position and updates match or
mismatch counts for many references at once. If a queried base is common at
that position, the matching list may contain most of the library, so it is
cheaper to update only the shorter mismatch list and treat all other references
as implicit matches. If the queried base is rare, the program uses the shorter
match list to update the matched-base counts for the corresponding references.
This adaptive choice matters because most
positions are conserved, while a small number of variant positions distinguish
the references.

For example, consider a toy 5-bp read (`AAAAA`) used only to illustrate the counting
logic. Suppose the library contains 1,000 reference sequences and this 5-bp read
is tested at start position 50, so its five bases are compared with positions
50-54 in all 1,000 references at the same time. At positions 50, 51, and 54,
the query bases are common in the library, so the program follows the mismatch
rule: it updates the mismatched-base counts for references that do not have the
query base, while the remaining references receive implicit matched bases. At
positions 52 and 53, the query bases are rare, so the program follows the match
rule: it updates the matched-base counts for references that do have the query
base, while references absent from that match list are treated as mismatches
when the final count is calculated.

![TEST SVG](data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxMDAiIGhlaWdodD0iMTAwIj48Y2lyY2xlIGN4PSI1MCIgY3k9IjUwIiByPSI0MCIgc3Ryb2tlPSJibGFjayIgc3Ryb2tlLXdpZHRoPSIzIiBmaWxsPSJyZWQiIC8+PC9zdmc+)

```text
index data at position 50:
  'matched' => 50 => 'A' => ['seq 001', 'seq 002', ..., 'seq 137', ... ]   (around 1000 items)
  'matched' => 50 => 'T' => ['seq 101', 'seq 302', ... ]   (very few items)
  'matched' => 50 => 'G' => ['seq 221', 'seq 252', ... ]   (very few items)
  'matched' => 50 => 'C' => ['seq 231', 'seq 288', ... ]   (very few items)

  'mismatched' => 50 => 'A' => ['seq 101', 'seq 302', ... ]   (very few items)
  'mismatched' => 50 => 'T' => ['seq 001', 'seq 002', ..., 'seq 137', ... ]   (around 1000 items)
  'mismatched' => 50 => 'G' => ['seq 001', 'seq 002', ..., 'seq 137', ... ]   (around 1000 items)
  'mismatched' => 50 => 'C' => ['seq 001', 'seq 002', ..., 'seq 137', ... ]   (around 1000 items)
```

The batched update produces separate counters for every reference. To see how
one counter is interpreted, consider reference 137. At the three mismatch-rule
positions, suppose reference 137 mismatches the read at position 51 only. The
program therefore stores one explicit mismatch for reference 137, and the other
two mismatch-rule positions are counted as implicit matches. At the two
match-rule positions, suppose reference 137 appears in the match list at
position 52 but not at position 53. The program therefore stores one explicit
match from the match-rule positions. 

```text
using index at position 50 for mapping the top 5-bp read

```

```text
numbers of matched and mismatched bases of this 5bp read:
  seq 001: 2 matched out of 2, 0 mismatched out of 3
  seq 002: 1 matched out of 2, 0 mismatched out of 3
  ...
  seq 137: 1 matched out of 2, 1 mismatched out of 3
  ...
```

The total matched count for reference 137
is:

```text
explicit matches from match-rule positions
+ implicit matches from mismatch-rule positions
= 1 + (3 - 1)
= 3 matched bases
```

The read overlaps five reference bases in this example, so the total mismatch
count for reference 137 is:

```text
aligned bases - matched bases = 5 - 3 = 2 mismatched bases
```

These two mismatches correspond to position 51, found through the mismatch
list, and position 53, inferred because the reference did not appear in the
rare-base match list. Thus, both rules contribute to the same final
match/mismatch totals while avoiding updates to long posting lists.

For each candidate start, the program accumulates two scores for every active
reference: the number of aligned bases that match the read and the number that
mismatch it. After all useful bases for that start have been considered, the
program compares the candidate score with the current best score for that same
reference. If the candidate has more matched bases, it becomes the new best
mapping location for that reference. If the matched-base count is tied, the
program uses the smaller mismatch count as the tie-breaker. In this way, each
reference sequence retains its own best start position, matched-base count, and
mismatched-base count for the read.

# Early Drop-out of Impossible Best Alignements
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
