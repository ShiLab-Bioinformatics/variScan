package util

import "sort"
import "fmt"

type SortedTable struct{
	SortedEntries []int
}

func CreateSortedTable() *SortedTable{
	ret := new(SortedTable)
	ret.SortedEntries = make([]int,0)

	return ret
}

func (me * SortedTable)Append(pos int){
	me.SortedEntries = append(me.SortedEntries, pos)
}

func (me * SortedTable)Sort(){
	sort.Sort(sort.IntSlice(me.SortedEntries))
}

func (me * SortedTable)LessOrEqual(pos int) int{
	var idxmin, idxmax int
	idxmin = 0
	idxmax = len(me.SortedEntries)

	if len(me.SortedEntries)<1 || me.SortedEntries[0] > pos{
		return -1
	}

	startidx:=0
	for {
		idxmid :=idxmin + (idxmax - idxmin) / 2
		if idxmid > idxmax{
			idxmid = idxmax
		}
		if me.SortedEntries[idxmid] < pos{
			if idxmax==idxmin+1{
				startidx = idxmin
				if idxmax < len(me.SortedEntries) && me.SortedEntries[idxmax] <= pos{ startidx = idxmax }
				break
			}
			idxmin = idxmid
		}else if me.SortedEntries[idxmid] > pos{
			idxmax = idxmid - 1
		}else{
			startidx = idxmid
			break
		}
		//fmt.Printf("%d ~ %d ~ %d\n", idxmin, idxmid, idxmax)

		if idxmin >= idxmax{
			startidx = idxmin
			break // now idxmid m
		}
	}

	return me.SortedEntries[startidx]
}

func (me * SortedTable)Range(pos_start, pos_stop int) []int{
	var idxmin, idxmax int

	idxmin = 0
	idxmax = len(me.SortedEntries)

	if me.SortedEntries[0] > pos_stop{
		return make([]int, 0)
	}

	if me.SortedEntries[len(me.SortedEntries) - 1] < pos_start{
		return make([]int, 0)
	}



	startidx:=0
	stopidx :=0
	for {
		idxmid :=idxmin + (idxmax - idxmin) / 2
		if idxmid > idxmax{
			idxmid = idxmax
		}
		if me.SortedEntries[idxmid] < pos_start{
			if idxmax==idxmin+1{
				startidx = idxmin
				break
			}
			idxmin = idxmid
		}else if me.SortedEntries[idxmid] > pos_start{
			idxmax = idxmid - 1
		}else{
			startidx = idxmid
			break
		}
		//fmt.Printf("%d ~ %d ~ %d\n", idxmin, idxmid, idxmax)

		if idxmin >= idxmax{
			startidx = idxmin
			break // now idxmid m
		}
	}


	idxmin = 0
	idxmax = len(me.SortedEntries)
	for {
		idxmid :=idxmin + (idxmax - idxmin) / 2

		if idxmid > idxmax{
			idxmid = idxmax
		}

		if me.SortedEntries[idxmid] < pos_stop{
			idxmin = idxmid + 1
		}else if me.SortedEntries[idxmid] > pos_stop{
			if idxmax==idxmin+1{
				stopidx = idxmax
			}
			idxmax = idxmid
		}else{
			stopidx = idxmid
			break
		}

		if idxmin >= idxmax{
			stopidx = idxmin
			break // now idxmid m
		}
	}

	if stopidx>=len(me.SortedEntries){
		stopidx = len(me.SortedEntries)-1
	}

	if stopidx>=startidx{
		ret := make([]int,0)
		for _, ri := range me.SortedEntries[startidx:stopidx+1]{
			if ri >= pos_start && ri<=pos_stop{
				ret = append(ret, ri)
			}
		}
		return ret
	}
	return make([]int, 0)
}

func (me * SortedTable)PrintEntries(){
	fmt.Printf("\n\n=========== Entries = %d ===========\n", len(me.SortedEntries))
	rlen := 0
	for _, ee := range me.SortedEntries{
		fmt.Printf("% 7d ", ee)
		rlen ++
		if rlen > 13{
			rlen = 0
			fmt.Println("")
		}
	}
}
