package util

import "sort"

func IfElseInt(con bool, i1, i2 int) int{
	if(con){
		return i1
	}
	return i2
}

func IfElseStr(con bool, i1, i2 string) string{
	if(con){
		return i1
	}
	return i2
}

func AbsInt(i int) int{
	if i>=0{
		return i
	}
	return -i
}

func MinInt(i,j int) int{
	if i>j{return j}
	return i
}


func MaxFloat64(i,j float64) float64{
	if i<j{return j}
	return i
}


func MaxInt(i,j int) int{
	if i<j{return j}
	return i
}


func IntersectInt(s1, s2 []int) []int{
	ret := make([]int,0)
	tab1 := make(map[int] bool)
	for _, r1 := range s1{
		tab1[r1]=true
	}

	for _, r2 := range s2{
		if tab1[r2]{
			ret = append(ret, r2)
		}
	}

	return ret
}

func IntersectStr(s1, s2 []string) []string{
	ret := make([]string,0)
	tab1 := make(map[string] bool)
	for _, r1 := range s1{
		tab1[r1]=true
	}

	for _, r2 := range s2{
		if tab1[r2]{
			ret = append(ret, r2)
		}
	}

	return ret
}

func UnionInt(s1, s2 []int) []int{
	tab :=  make(map[int] bool)
	ret := make([]int,0)

	for _, r := range s1{
		tab[r]=true
	}

	for _, r := range s2{
		tab[r]=true
	}

	for r := range tab{
		ret=append(ret,r)
	}
	return ret
}

func UnionStr(s1, s2 []string) []string{
	tab :=  make(map[string] bool)
	ret := make([]string,0)

	for _, r := range s1{
		tab[r]=true
	}

	for _, r := range s2{
		tab[r]=true
	}

	for r := range tab{
		ret=append(ret,r)
	}
	return ret
}
type AsGap [][2]int

func (a AsGap) Len() int           { return len(a) }
func (a AsGap) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a AsGap) Less(i, j int) bool { return a[i][0] < a[j][0] }

// Gaps are open gaps. Close gap of [1, 99] is [1, 100) in input.
// Output is also open gaps.
// gapA and gapB MUST be sorted by coordinates and non-overlapping.
func OverlapGaps(gapA [][2]int, gapB [][2]int) [][2]int{
	ret := make([][2]int,0)
	iA:=0
	iB:=0

	poscur := MaxInt(gapA[0][0], gapB[0][0])
	for{
		cur_end := MinInt(gapA[iA][1], gapB[iB][1])
		if cur_end > poscur { ret = append(ret,[2]int{poscur, cur_end}) }
		// can be both:
		if cur_end == gapB[iB][1] {iB++}
		if cur_end == gapA[iA][1] {iA++}
		if iB >=len(gapB) || iA >=len(gapA){break}

		poscur = MaxInt(gapA[iA][0], gapB[iB][0])
	}
	return ret
}

func MergedGaps(gaps [][2]int) [][2]int{
	ret := make([][2]int,0)
	sort.Sort( AsGap(gaps) )

	ret = append( ret, gaps[0] )
	for _, g1 := range gaps[1:]{
		if ret[ len(ret) -1 ][1] >= g1[1]{
			continue
		}else if g1[0] > ret[ len(ret) -1 ][1]{
			ret = append(ret, g1)
		}else{
			ret[ len(ret) -1 ][1] = g1[1]
		}
	}
	return ret
}

func UniqueStr(strs[]string) []string{
	smap := make(map[string]bool)
	for _,s := range strs{
		smap[s]=true
	}

	ret := make([]string, 0)
	for s,_ := range smap{
		ret=append(ret, s)
	}

	return ret
}

func UniqueInt(strs[]int) []int{
	smap := make(map[int]bool)
	for _,s := range strs{
		smap[s]=true
	}

	ret := make([]int, 0)
	for s,_ := range smap{
		ret=append(ret, s)
	}

	return ret
}

