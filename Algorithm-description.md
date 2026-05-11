# Positional Index Alignment Algorithm

This program aligns short reads against a library of thousands of closely
related reference sequences. 

The algorithm does not perform full dynamic programming alignment and does not
model insertions or deletions. Instead, it searches for the best ungapped
placement of each read on each reference. For every reference sequence, the
best placement is the start position with the highest number of matching bases.
At that placement, the program reports:

- the start position,
- the number of matched bases, and
- the number of mismatched bases.

Ties in matched-base count are resolved by choosing the placement with fewer
mismatches. This alignment model is especially suitable for reference libraries
where sequences are similar in length and differ mainly by substitutions.

## Core Idea

The key observation is that the reference sequences are highly similar and are
usually the same length. A naive implementation would compare every read with
every reference at every possible start position. For `R` references, read
length `L`, and `S` candidate starts, this requires roughly:

```text
R * L * S
```

base comparisons per read. With thousands of near-identical references, much of
this work is redundant because the same conserved bases are checked repeatedly.

This algorithm changes the unit of work. Instead of asking:

```text
How well does this read match reference 1?
How well does this read match reference 2?
How well does this read match reference 3?
...
```

it asks:

```text
At this aligned position, which references contain the query base?
```

One lookup can therefore update the score of many references at once.

## Positional Inverted Index

To support this batched scoring, the program builds a positional inverted index
over the equal-length reference sequences. For each reference position and each
base `A/C/G/T`, the index stores two posting lists:

- reference sequence IDs that match that base at the position, and
- reference sequence IDs that mismatch that base at the position.

Conceptually, the index has the following structure:

```text
[match or mismatch] -> [reference position] -> [base] -> [reference sequence IDs]
```

Suppose the library contains 1,000 reference sequences, and suppose that, at base-position 50, most reference sequences
contain `A`, while only a small number contain `T`, `G`, or `C`. The index for
that position would look conceptually like this:

```text
'matched'    => 'A' => ['seq 001', 'seq 002', ...]      (nearly 1000 items)
'matched'    => 'T' => ['seq 101', 'seq 302', ...]      (very few items)
'matched'    => 'G' => ['seq 221', 'seq 252', ...]      (very few items)
'matched'    => 'C' => ['seq 231', 'seq 288', ...]      (very few items)

'mismatched' => 'A' => ['seq 101', 'seq 221', 'seq 231', ...]        (very few items)
'mismatched' => 'T' => ['seq 001', 'seq 002', ...]      (nearly 1000 items)
'mismatched' => 'G' => ['seq 001', 'seq 002', ...]      (nearly 1000 items)
'mismatched' => 'C' => ['seq 001', 'seq 002', ...]      (nearly 1000 items)
```

The `matched A` list is long because most reference sequences have `A` at this position.
The `mismatched A` list is short because only the reference sequences with another base
at this position fail to match `A`. The reverse is true for rare bases such as
`T`, `G`, and `C`: their match lists are short, while their mismatch lists are
long.

This redundancy is the source of the speed-up. When one posting list is long,
the complementary list is short, and the algorithm can update only the shorter
one.

## Simultaneous Match and Mismatch Counting

During alignment, a read is placed at a __candidate start__ position. Only the
overlapping part of the read and reference is evaluated, so partial overlaps
near reference ends are also supported.

For each overlapping query base, the program looks up the corresponding
reference position in the index. It then chooses the shorter of two equivalent
ways to score the base:

- **Mismatch mode:** If the query base is common at this position, the match
  list is long and the mismatch list is short. The program updates only the
  references in the mismatch list. All other reference sequences are treated as implicit
  matches.
- **Match mode:** If the query base is rare at this position, the match list is
  short and the mismatch list is long. The program updates only the reference sequences
  in the match list. All other reference sequences are treated as implicit mismatches
  when the final score is calculated.

This adaptive choice is important because most positions are conserved, while a
small number of variant positions distinguish the reference sequences.

## Counting Example

Consider a toy 5-bp read:

```text
AAAAA
```

Suppose the library contains 1,000 reference sequences and this read is tested
at start position 50. The five query bases are therefore compared with reference
positions 50 to 54 across all 1,000 reference sequences at the same time.

Assume the query base `A` is common at positions 50, 51, and 54 in these 1000 reference sequences. For those 3
positions, the program uses `mismatch mode`: it records only the few references
that do not have `A`, and it treats all other reference sequences as implicit matches.

Now assume the query base `A` is rare at positions 52 and 53 in these 1000 reference sequences. For those 2
positions, the program uses `match mode`: it records only the few reference sequences that
do have `A`, and it treats all other references as implicit mismatches.

### Common Query Base

At position 50, most references have `A`. Therefore, for the first base of the
read, the `matched A` list is much longer than the `mismatched A` list:

```text
'matched'    => 'A' => ['seq 001', 'seq 002', ...]      (nearly 1000 items)
'mismatched' => 'A' => ['seq 101', 'seq 221', ...]      (very few items)
```

The algorithm uses mismatch mode and updates only the few reference sequences in
`mismatched A`. Every reference sequence absent from that short list receives an implicit
match for this read base.

### Rare Query Base

At position 52, suppose most reference sequences have `C`, while only a few reference sequences
have `A`. The relevant index entries are therefore:

```text
'matched'    => 'A' => ['seq 137', 'seq 261', ...]      (very few items)
'matched'    => 'C' => ['seq 001', 'seq 002', ...]      (nearly 1000 items)

'mismatched' => 'A' => ['seq 001', 'seq 002', ...]      (nearly 1000 items)
'mismatched' => 'C' => ['seq 137', 'seq 285', ...]      (very few items)
```

For the third base of the read, the query base is still `A`, but `A` is rare at
this reference sequence position. The `matched A` list is therefore much shorter than the
`mismatched A` list. The algorithm uses match mode and updates only the
reference sequences in `matched A`. Reference sequences absent from this short list are counted as
mismatches for this position when the final total is computed.

### Combining Explicit and Implicit Counts

The program combines explicit and implicit counts to obtain the final matched/mismatched counts of this read to each reference sequences at that candidate start. As described above,

```text
2 in the 5 read bases were compared with the reference sequences on the match mode.
3 in the 5 read bases were compared with the reference sequences on the mismatch mode.
```

For example, suppose `seq 137` has:

```text
1 explicit match in the 2 match-mode positions.
1 explicit mismatch in the 3 mismatch-mode positions.
```

These explicit counts also determine the implicit counts:

```text
match-mode positions:
  2 total positions - 1 explicit match = 1 implicit mismatch

mismatch-mode positions:
  3 total positions - 1 explicit mismatch = 2 implicit matches
```

Hence, the total matched-base count in `seq 137` is:

```text
explicit matches from match-mode positions
+ implicit matches from mismatch-mode positions
= 1 + 2
= 3 matched bases

```

The total mismatched-base count is:

```text
explicit mismatches from mismatch-mode positions
+ implicit mismatches from match-mode positions
= 1 + 1
= 2 mismatched bases
```

Both modes therefore contribute to the same final match and mismatch totals,
while avoiding updates to long posting lists.

## Choosing the Best Start for Each Reference

For each candidate start, the program accumulates a score for every 
reference sequence. After all overlapping bases for that candidate start have been compared between read and all reference sequences, each score is compared with the best scores for each reference sequence.

A candidate start becomes the new best placement for a reference if:

- it has more matched bases than the current best placement, or
- it has the same number of matched bases but fewer mismatched bases.

In this way, every reference sequence keeps its own best start position,
matched-base count, and mismatched-base count for the read.

## Shorter Reference Sequences

Reference sequences shorter than the longest sequence (typically the majority) are processed separately using direct scanning, because they cannot share the same fixed-position index.

For these shorter references, the program evaluates every possible ungapped alignment start position, directly counts matched and mismatched bases, and retains the placement with the highest number of matched bases. The same tie-breaking rule is applied: if two placements have the same number of matched bases, the placement with fewer mismatches is preferred.

Even when the number of shorter reference sequences exceeds the number of longest references, the alignment results remain correct, while the algorithm becomes slower.

## Parallel Execution and Output

The program is parallelized across reads using worker goroutines. The reference
index is built once, shared read-only by all workers, and then reused for
independent read alignments. This is effective because reads do not depend on
one another after the reference index has been constructed.
