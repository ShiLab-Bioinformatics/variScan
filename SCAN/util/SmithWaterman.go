package util

import (
	"strings"
	"math"
	"bytes"
	"strconv"
	"fmt"
)

// AlignmentResult holds the results of the Smith-Waterman alignment.
type AlignmentResult struct {
	Score      int
	SequenceA  string
	SequenceB  string
	StartA     int // Start index in the original Sequence A
	StartB     int // Start index in the original Sequence B
}

const (
	STOP = 0
	DIAG = 1 // Diagonal move (match/mismatch)
	UP   = 2 // Up move (gap in sequence B)
	LEFT = 3 // Left move (gap in sequence A)
)



// SmithWaterman performs local sequence alignment using a traceback matrix for correctness.
func SmithWaterman(seqA, seqB string, matchScore, mismatchScore, gapCreationPenalty, gapExtensionPenalty int) AlignmentResult {
	lenA := len(seqA)
	lenB := len(seqB)

	// 1. Initialize the scoring and traceback matrices.
	H := make([][]int, lenA+1)
	traceback := make([][]int, lenA+1)
	for i := range H {
		H[i] = make([]int, lenB+1)
		traceback[i] = make([]int, lenB+1) // Defaults to STOP (0)
	}

	maxScore := 0
	maxI, maxJ := 0, 0

	// 2. Fill the matrices.
	for i := 1; i <= lenA; i++ {
		for j := 1; j <= lenB; j++ {
			// Calculate the score for a diagonal move.
			diagScore := H[i-1][j-1]
			if seqA[i-1] == seqB[j-1] {
				diagScore += matchScore
			} else {
				diagScore += mismatchScore
			}

			// Calculate scores for moves from up and left.
			upIsExtension := traceback[i-1][j] != DIAG
			upScore := H[i-1][j] + IfElseInt(upIsExtension,gapExtensionPenalty,gapCreationPenalty)

			leftIsExtension := traceback[i][j-1] != DIAG
			leftScore := H[i][j-1] + IfElseInt(leftIsExtension,gapExtensionPenalty,gapCreationPenalty)

			// Determine the maximum score and record the direction.
			// The order of checks (diag, up, left) sets the priority if scores are equal.
			currentMax := 0
			direction := STOP

			if diagScore > currentMax {
				currentMax = diagScore
				direction = DIAG
			}
			if upScore > currentMax {
				currentMax = upScore
				direction = UP
			}
			if leftScore > currentMax {
				currentMax = leftScore
				direction = LEFT
			}

			H[i][j] = currentMax
			traceback[i][j] = direction

			// Update the overall maximum score found so far.
			if currentMax > maxScore {
				maxScore = currentMax
				maxI, maxJ = i, j
			}
		}
	}

	// If no positive alignment was found, return an empty result.
	if maxScore == 0 {
		return AlignmentResult{Score: 0, SequenceA: "", SequenceB: ""}
	}

	// 3. Perform traceback using the direction matrix.
	alignedA := strings.Builder{}
	alignedB := strings.Builder{}
	i, j := maxI, maxJ

	for traceback[i][j] != STOP {
		switch traceback[i][j] {
		case DIAG:
			alignedA.WriteByte(seqA[i-1])
			alignedB.WriteByte(seqB[j-1])
			i--
			j--
		case UP:
			alignedA.WriteByte(seqA[i-1])
			alignedB.WriteByte('-')
			i--
		case LEFT:
			alignedA.WriteByte('-')
			alignedB.WriteByte(seqB[j-1])
			j--
		}
	}

	return AlignmentResult{
		Score:      maxScore,
		SequenceA:  reverseString(alignedA.String()),
		SequenceB:  reverseString(alignedB.String()),
		StartA:     i, // Final i and j are the start indices.
		StartB:     j,
	}
}

// Helper function to reverse a string.
func reverseString(s string) string {
	runes := []rune(s)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}









// SmithWaterman performs local sequence alignment using a traceback matrix for correctness.
func SmithWaterman3P(seqA, seqB string, matchScore, mismatchScore, gapPenalty int) AlignmentResult {
	lenA := len(seqA)
	lenB := len(seqB)

	// 1. Initialize the scoring and traceback matrices.
	H := make([][]int, lenA+1)
	traceback := make([][]int, lenA+1)
	for i := range H {
		H[i] = make([]int, lenB+1)
		traceback[i] = make([]int, lenB+1) // Defaults to STOP (0)
	}

	maxScore := 0
	maxI, maxJ := 0, 0

	// 2. Fill the matrices.
	for i := 1; i <= lenA; i++ {
		for j := 1; j <= lenB; j++ {
			// Calculate the score for a diagonal move.
			diagScore := H[i-1][j-1]
			if seqA[i-1] == seqB[j-1] {
				diagScore += matchScore
			} else {
				diagScore += mismatchScore
			}

			// Calculate scores for moves from up and left.
			upScore := H[i-1][j] + gapPenalty
			leftScore := H[i][j-1] + gapPenalty

			// Determine the maximum score and record the direction.
			// The order of checks (diag, up, left) sets the priority if scores are equal.
			currentMax := 0
			direction := STOP

			if diagScore > currentMax {
				currentMax = diagScore
				direction = DIAG
			}
			if upScore > currentMax {
				currentMax = upScore
				direction = UP
			}
			if leftScore > currentMax {
				currentMax = leftScore
				direction = LEFT
			}

			H[i][j] = currentMax
			traceback[i][j] = direction

			// Update the overall maximum score found so far.
			if currentMax > maxScore {
				maxScore = currentMax
				maxI, maxJ = i, j
			}
		}
	}

	// If no positive alignment was found, return an empty result.
	if maxScore == 0 {
		return AlignmentResult{Score: 0, SequenceA: "", SequenceB: ""}
	}

	// 3. Perform traceback using the direction matrix.
	alignedA := strings.Builder{}
	alignedB := strings.Builder{}
	i, j := maxI, maxJ

	for traceback[i][j] != STOP {
		switch traceback[i][j] {
		case DIAG:
			alignedA.WriteByte(seqA[i-1])
			alignedB.WriteByte(seqB[j-1])
			i--
			j--
		case UP:
			alignedA.WriteByte(seqA[i-1])
			alignedB.WriteByte('-')
			i--
		case LEFT:
			alignedA.WriteByte('-')
			alignedB.WriteByte(seqB[j-1])
			j--
		}
	}

	return AlignmentResult{
		Score:      maxScore,
		SequenceA:  reverseString(alignedA.String()),
		SequenceB:  reverseString(alignedB.String()),
		StartA:     i, // Final i and j are the start indices.
		StartB:     j,
	}
}



type ScoringParams struct {
	MatchScore      int
	MismatchPenalty int
	GapCreatePenalty int // Cost to open a new gap
	GapExtPenalty   int // Cost to extend an existing gap
}


const negInf = math.MinInt32 / 2 // Use half min int to prevent overflow during additions

// max3 returns the max of three ints.
func max3(a, b, c int) int {
	if a >= b && a >= c {
		return a
	}
	if b >= a && b >= c {
		return b
	}
	return c
}

// max2 returns the max of two ints.
func max2(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// AffineGlobalAlignment performs global alignment (End-to-End) so no base is left unmapped.
func AffineGlobalAlignment(query, reference string, params ScoringParams) (cigar string, score int) {
	r1 := []rune(query)
	r2 := []rune(reference)
	n := len(r1)
	m := len(r2)

	// 2. Initialize Matrices
	// M: Optimal score ending in a Match/Mismatch
	// X: Optimal score ending with a gap in Seq2 (Insertion in Seq1)
	// Y: Optimal score ending with a gap in Seq1 (Insertion in Seq2)
	M := make([][]int, n+1)
	X := make([][]int, n+1)
	Y := make([][]int, n+1)

	for i := range M {
		M[i] = make([]int, m+1)
		X[i] = make([]int, m+1)
		Y[i] = make([]int, m+1)
	}

	// 3. Boundary Initialization (Crucial for Global Alignment)
	// Unlike Smith-Waterman (which zeros these), we must accumulate gap penalties 
	// to force the alignment to start at index 0.

	M[0][0] = 0
	X[0][0] = negInf
	Y[0][0] = negInf

	for i := 1; i <= n; i++ {
		// Cost to have a gap of length i at the start of query
		cost := -(params.GapCreatePenalty + (i * params.GapExtPenalty))
		M[i][0] = negInf // Cannot match against nothing
		X[i][0] = cost   // Gap in reference (vertical move)
		Y[i][0] = negInf // Invalid state
	}

	for j := 1; j <= m; j++ {
		// Cost to have a gap of length j at the start of reference
		cost := -(params.GapCreatePenalty + (j * params.GapExtPenalty))
		M[0][j] = negInf
		X[0][j] = negInf
		Y[0][j] = cost   // Gap in query (horizontal move)
	}

	// 4. Fill Matrices
	for i := 1; i <= n; i++ {
		for j := 1; j <= m; j++ {
			// --- Calculate X (Vertical / Insertion in Seq1) ---
			// Extend existing X gap OR Start new gap from M
			// Note: We assume 'Create' is paid once per gap group.
			extendX := X[i-1][j] - params.GapExtPenalty
			openX := M[i-1][j] - params.GapCreatePenalty - params.GapExtPenalty
			X[i][j] = max2(extendX, openX)

			// --- Calculate Y (Horizontal / Insertion in Seq2) ---
			extendY := Y[i][j-1] - params.GapExtPenalty
			openY := M[i][j-1] - params.GapCreatePenalty - params.GapExtPenalty
			Y[i][j] = max2(extendY, openY)

			// --- Calculate M (Match or Mismatch) ---
			score := -params.MismatchPenalty
			if r1[i-1] == r2[j-1] {
				score = params.MatchScore
			}

			// Can come from previous Match, Close X, or Close Y
			fromM := M[i-1][j-1] + score
			fromX := X[i-1][j-1] + score
			fromY := Y[i-1][j-1] + score

			M[i][j] = max3(fromM, fromX, fromY)
		}
	}

	// 5. Traceback
	// Start from the very bottom-right corner (Global alignment requirement)
	i, j := n, m

	// Determine which matrix produced the final best score
	finalScore := max3(M[n][m], X[n][m], Y[n][m])

	// State tracker: 0=M, 1=X, 2=Y
	state := 0
	if finalScore == X[n][m] {
		state = 1
	} else if finalScore == Y[n][m] {
		state = 2
	}

	var align1, align2 strings.Builder

	for i > 0 || j > 0 {
		if state == 0 { // In M (Match/Mismatch) state
			// Determine where we came from
			score := -params.MismatchPenalty
			if r1[i-1] == r2[j-1] {
				score = params.MatchScore
			}

			// Recalculate sources to find path
			// Check if we came from M
			if i > 0 && j > 0 && M[i][j] == M[i-1][j-1]+score {
				align1.WriteRune(r1[i-1])
				align2.WriteRune(r2[j-1])
				i--
				j--
				state = 0
			} else if i > 0 && j > 0 && M[i][j] == X[i-1][j-1]+score {
				align1.WriteRune(r1[i-1])
				align2.WriteRune(r2[j-1])
				i--
				j--
				state = 1 // Moved to X
			} else {
				align1.WriteRune(r1[i-1])
				align2.WriteRune(r2[j-1])
				i--
				j--
				state = 2 // Moved to Y
			}

		} else if state == 1 { // In X (Vertical gap) state
			align1.WriteRune(r1[i-1])
			align2.WriteRune('-')
			// Did we extend X or open X from M?
			extendCost := X[i-1][j] - params.GapExtPenalty
			if X[i][j] == extendCost {
				state = 1 // Stay in X
			} else {
				state = 0 // Close gap, go to M
			}
			i--

		} else if state == 2 { // In Y (Horizontal gap) state
			align1.WriteRune('-')
			align2.WriteRune(r2[j-1])

			// Did we extend Y or open Y from M?
			extendCost := Y[i][j-1] - params.GapExtPenalty
			if Y[i][j] == extendCost {
				state = 2 // Stay in Y
			} else {
				state = 0 // Close gap, go to M
			}
			j--
		}
	}

        score=finalScore

	AlignSeqQuery:= reverseString(align1.String())
	AlignSeqRef  := reverseString(align2.String())
	if len(AlignSeqQuery) != len(AlignSeqRef) {panic("Internal error: unequal cigar lengths")}
	cigar = ""
	if false{
		println(AlignSeqQuery)
		println(AlignSeqRef)
	}
	tmpi := 0
	last_nch := 'x'
	for Qi, Qch := range AlignSeqQuery{
		nch := 'x'
		Rch := AlignSeqRef[Qi]
		if Qch == '-' && Rch == '-' {panic("BOTH UNMATCHED")}
		if Qch == '-' {
			nch = 'D'
		}else if Rch == '-' {
			nch = 'I'
		}else{
			nch = 'M'
		}
		if last_nch != nch{
			if last_nch!='x'{cigar += fmt.Sprintf("%d%c",tmpi, last_nch)}
			tmpi = 0
			last_nch = nch
		}
		tmpi ++
	}
	if last_nch!='x' {cigar += fmt.Sprintf("%d%c",tmpi, last_nch)}
	return
}

// Scores and Penalties are all positive.
func GlobalAlignmentCIGAR(query_seq, reference_seq string, matchScore, mismatchPenalty, gapOpenPenalty, gapExtendPenalty int) (string, int) {
	seq2:= query_seq
	seq1:= reference_seq
	n := len(seq1)
	m := len(seq2)

	const negInf = -1 << 60

	// DP matrices: M = match/mismatch, Ix = gap in seq2 (D), Iy = gap in seq1 (I)
	M := make([][]int, m+1)
	Ix := make([][]int, m+1)
	Iy := make([][]int, m+1)
	for i := 0; i <= m; i++ {
		M[i] = make([]int, n+1)
		Ix[i] = make([]int, n+1)
		Iy[i] = make([]int, n+1)
	}

	type state uint8
	const (
		stateM state = iota
		stateIx
		stateIy
	)

	// traceback state matrices
	prevStateM := make([][]state, m+1)
	prevStateIx := make([][]state, m+1)
	prevStateIy := make([][]state, m+1)
	for i := 0; i <= m; i++ {
		prevStateM[i] = make([]state, n+1)
		prevStateIx[i] = make([]state, n+1)
		prevStateIy[i] = make([]state, n+1)
	}

	// Initialize matrices
	for i := 0; i <= m; i++ {
		for j := 0; j <= n; j++ {
			M[i][j] = negInf
			Ix[i][j] = negInf
			Iy[i][j] = negInf
		}
	}

	M[0][0] = 0
	prevStateM[0][0] = stateM

	// First row: gaps in seq2 => Ix (D)
	for j := 1; j <= n; j++ {
		if j == 1 {
			Ix[0][j] = M[0][0] - (gapOpenPenalty + gapExtendPenalty)
			prevStateIx[0][j] = stateM
		} else {
			Ix[0][j] = Ix[0][j-1] - gapExtendPenalty
			prevStateIx[0][j] = stateIx
		}
	}

	// First column: gaps in seq1 => Iy (I)
	for i := 1; i <= m; i++ {
		if i == 1 {
			Iy[i][0] = M[0][0] - (gapOpenPenalty + gapExtendPenalty)
			prevStateIy[i][0] = stateM
		} else {
			Iy[i][0] = Iy[i-1][0] - gapExtendPenalty
			prevStateIy[i][0] = stateIy
		}
	}

	// Fill DP matrices
	for i := 1; i <= m; i++ {
		for j := 1; j <= n; j++ {
			// score for aligning seq2[i-1] with seq1[j-1]
			var pairScore int
			if seq2[i-1] == seq1[j-1] {
				pairScore = matchScore
			} else {
				pairScore = -mismatchPenalty
			}

			// M[i][j]
			best := M[i-1][j-1]
			prev := stateM
			if Ix[i-1][j-1] > best {
				best = Ix[i-1][j-1]
				prev = stateIx
			}
			if Iy[i-1][j-1] > best {
				best = Iy[i-1][j-1]
				prev = stateIy
			}
			M[i][j] = best + pairScore
			prevStateM[i][j] = prev

			// Ix[i][j] : gap in seq2 (D)
			open := M[i][j-1] - (gapOpenPenalty + gapExtendPenalty)
			extend := Ix[i][j-1] - gapExtendPenalty
			if open > extend {
				Ix[i][j] = open
				prevStateIx[i][j] = stateM
			} else {
				Ix[i][j] = extend
				prevStateIx[i][j] = stateIx
			}

			// Iy[i][j] : gap in seq1 (I)
			open = M[i-1][j] - (gapOpenPenalty + gapExtendPenalty)
			extend = Iy[i-1][j] - gapExtendPenalty
			if open > extend {
				Iy[i][j] = open
				prevStateIy[i][j] = stateM
			} else {
				Iy[i][j] = extend
				prevStateIy[i][j] = stateIy
			}
		}
	}

	// Choose best end state
	endI, endJ := m, n
	endState := stateM
	bestScore := M[m][n]
	if Ix[m][n] > bestScore {
		bestScore = Ix[m][n]
		endState = stateIx
	}
	if Iy[m][n] > bestScore {
		bestScore = Iy[m][n]
		endState = stateIy
	}

	// Traceback to build CIGAR (in reverse)
	ops := make([]byte, 0, m+n)
	i, j, st := endI, endJ, endState
	for i > 0 || j > 0 {
		switch st {
		case stateM:
			// Diagonal move: match/mismatch
			ops = append(ops, 'M')
			ps := prevStateM[i][j]
			i--
			j--
			st = ps
		case stateIx:
			// Left move: deletion in query => 'D'
			ops = append(ops, 'D')
			ps := prevStateIx[i][j]
			j--
			st = ps
		case stateIy:
			// Up move: insertion in query => 'I'
			ops = append(ops, 'I')
			ps := prevStateIy[i][j]
			i--
			st = ps
		}
	}

	// Reverse ops
	for l, r := 0, len(ops)-1; l < r; l, r = l+1, r-1 {
		ops[l], ops[r] = ops[r], ops[l]
	}

	// Compress to CIGAR string
	var sb strings.Builder
	if len(ops) > 0 {
		cur := ops[0]
		count := 1
		for k := 1; k < len(ops); k++ {
			if ops[k] == cur {
				count++
			} else {
				sb.WriteString(fmt.Sprintf("%d%c", count, cur))
				cur = ops[k]
				count = 1
			}
		}
		sb.WriteString(fmt.Sprintf("%d%c", count, cur))
	}

	return sb.String(), bestScore
}





// SmithParams holds the penalties and scores.
// Note: ALL positive values
type SmithParams struct {
	MatchScore      int
	MismatchPenalty int
	GapOpenPenalty  int
	GapExtendPenalty int
}

// SmithResult holds the output of the algorithm.
type SmithResult struct {
	CIGAR     string
	RefStart  int // 0-based index start in Reference
	RefEnd    int // 0-based index end in Reference (exclusive)
	Score     int
	QuerySeq  string
	RefSeq    string
}

// MatrixCell holds scores for the three affine states.
type MatrixCell struct {
	mScore int // Match/Mismatch state
	iScore int // Insertion state (Gap in Reference, Extra in Query)
	dScore int // Deletion state (Gap in Query, Extra in Reference)
}

// SmithWatermanSemiGlobal performs alignment where the entire query must be mapped,
// but it can map to any substring of the reference.
func SmithWatermanSemiGlobal(query, ref string, params SmithParams) SmithResult {
	qLen := len(query)
	rLen := len(ref)

	// Initialize DP matrices
	// We use (qLen+1) x (rLen+1) grid
	dp := make([][]MatrixCell, qLen+1)
	for i := range dp {
		dp[i] = make([]MatrixCell, rLen+1)
	}

	// Constants to represent negative infinity for impossible states
	const negInf = math.MinInt32 / 2

	// 1. Initialization
	
	// [0][0]
	dp[0][0] = MatrixCell{mScore: 0, iScore: negInf, dScore: negInf}

	// Initialize Top Row (Reference Axis)
	// In Semi-Global (Local Ref), starting anywhere in Ref is free (0).
	// However, we cannot start in an Insert/Delete state at the very border easily 
	// without opening a gap, but usually Top Row is just 0 for M state.
	for j := 1; j <= rLen; j++ {
		dp[0][j] = MatrixCell{
			mScore: 0,      // Free start
			iScore: negInf, // Cannot have insertion in query at start of query
			dScore: negInf, // Initial deletion logic handled in loop if needed, but usually 0 or penalty
		}
	}

	// Initialize Left Column (Query Axis)
	// In Semi-Global (Global Query), skipping the start of Query costs points.
	for i := 1; i <= qLen; i++ {
		// Cost to have 'i' insertions (gaps in ref) at the start
		cost := params.GapOpenPenalty + (i * params.GapExtendPenalty)
		dp[i][0] = MatrixCell{
			mScore: negInf, // Cannot match empty ref
			iScore: -cost,  // Force gap creation
			dScore: negInf,
		}
	}

	// 2. Matrix Filling (Forward Pass)
	for i := 1; i <= qLen; i++ {
		for j := 1; j <= rLen; j++ {
			// Calculate match/mismatch score
			charScore := -params.MismatchPenalty
			if query[i-1] == ref[j-1] {
				charScore = params.MatchScore
			}

			// --- Calculate I_SCORE (Insertion in Query / Gap in Ref) ---
			// Extend an existing insertion or open a new one from a match/deletion
			iFromM := dp[i-1][j].mScore - (params.GapOpenPenalty + params.GapExtendPenalty)
			iFromI := dp[i-1][j].iScore - params.GapExtendPenalty
			// Depending on implementation, you might allow D->I, but typically strictly forbidden or expensive.
			// We will stick to standard Affine: M->I or I->I
			dp[i][j].iScore = max(iFromM, iFromI)

			// --- Calculate D_SCORE (Deletion in Query / Gap in Ref) ---
			// Extend existing deletion or open new one
			dFromM := dp[i][j-1].mScore - (params.GapOpenPenalty + params.GapExtendPenalty)
			dFromD := dp[i][j-1].dScore - params.GapExtendPenalty
			dp[i][j].dScore = max(dFromM, dFromD)

			// --- Calculate M_SCORE (Match/Mismatch) ---
			// Can come from diagonal M, I, or D state
			mFromM := dp[i-1][j-1].mScore + charScore
			mFromI := dp[i-1][j-1].iScore + charScore
			mFromD := dp[i-1][j-1].dScore + charScore
			
			dp[i][j].mScore = max(mFromM, max(mFromI, mFromD))
		}
	}

	// 3. Find Max Score in the Last Row (End of Query)
	// Because Query must be fully mapped, we only look at i == qLen.
	// We look for the best column j to end at.
	maxScore := negInf
	maxJ := 0
	endState := 0 // 0:M, 1:I, 2:D

	for j := 1; j <= rLen; j++ {
		cell := dp[qLen][j]
		
		// Check Match State
		if cell.mScore >= maxScore {
			maxScore = cell.mScore
			maxJ = j
			endState = 0
		}
		// Check Insertion State (Sequence ends with an insertion)
		if cell.iScore >= maxScore {
			maxScore = cell.iScore
			maxJ = j
			endState = 1
		}
		// Check Deletion State
		if cell.dScore >= maxScore {
			maxScore = cell.dScore
			maxJ = j
			endState = 2
		}
	}

	// 4. Traceback
	var cigarOps []string
	currI, currJ := qLen, maxJ
	currState := endState // Start tracking from the state that gave the max score

	// Note: We stop when currI is 0. 
	// currJ might be > 0 (meaning we didn't start at Ref[0], which is fine).
	for currI > 0 {
		// Based on current state, decide where we came from
		
		// If current is Match/Mismatch
		if currState == 0 {
			cigarOps = append(cigarOps, "M")
			
			// Recalculate where M came from
			score := dp[currI][currJ].mScore
			charScore := -params.MismatchPenalty
			if query[currI-1] == ref[currJ-1] {
				charScore = params.MatchScore
			}
			
			prevScore := score - charScore
			
			// Check which previous cell + charScore equals current score
			if prevScore == dp[currI-1][currJ-1].mScore {
				currState = 0
			} else if prevScore == dp[currI-1][currJ-1].iScore {
				currState = 1
			} else {
				currState = 2
			}
			currI--
			currJ--

		} else if currState == 1 { 
			// Current is Insertion (Gap in Ref, consumes Query)
			cigarOps = append(cigarOps, "I")
			
			score := dp[currI][currJ].iScore
			extendCost := params.GapExtendPenalty
			openCost := params.GapOpenPenalty
			
			// Did we extend or open?
			if score == dp[currI-1][currJ].iScore - extendCost {
				currState = 1
			} else if score == dp[currI-1][currJ].mScore - (openCost + extendCost) {
				currState = 0
			}
			currI-- // Consumed Query, Ref stays

		} else if currState == 2 {
			// Current is Deletion (Gap in Query, consumes Ref)
			cigarOps = append(cigarOps, "D")
			
			score := dp[currI][currJ].dScore
			extendCost := params.GapExtendPenalty
			openCost := params.GapOpenPenalty
			
			if score == dp[currI][currJ-1].dScore - extendCost {
				currState = 2
			} else if score == dp[currI][currJ-1].mScore - (openCost + extendCost) {
				currState = 0
			}
			currJ-- // Consumed Ref, Query stays
		}
	}

	// Reverse CIGAR ops
	for i, j := 0, len(cigarOps)-1; i < j; i, j = i+1, j-1 {
		cigarOps[i], cigarOps[j] = cigarOps[j], cigarOps[i]
	}

	// Compress CIGAR (e.g., M, M, M -> 3M)
	cigarStr := compressCigar(cigarOps)

	return SmithResult{
		CIGAR:    cigarStr,
		RefStart: currJ, // Where we stopped moving back in J is the start
		RefEnd:   maxJ,
		Score:    maxScore,
		QuerySeq: query,
		RefSeq:   ref,
	}
}

func compressCigar(ops []string) string {
	if len(ops) == 0 {
		return ""
	}
	var buf bytes.Buffer
	count := 1
	prev := ops[0]
	
	for i := 1; i < len(ops); i++ {
		if ops[i] == prev {
			count++
		} else {
			buf.WriteString(strconv.Itoa(count))
			buf.WriteString(prev)
			count = 1
			prev = ops[i]
		}
	}
	buf.WriteString(strconv.Itoa(count))
	buf.WriteString(prev)
	return buf.String()
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
