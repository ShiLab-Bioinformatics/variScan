package util

import "fmt"

type BucketedTable struct{
	dict map[string] map[int] []string
	granularity int
}


// granularity must be longer than the longest read length
func CreateBucketedTable(granularity int) * BucketedTable{
	ret := new(BucketedTable)
	ret.granularity = granularity
	ret.dict = make(map[string] map[int] []string)
	return ret
}

func (who * BucketedTable) Append(chroname string, pos int, value string) {
	ky1 := fmt.Sprintf("%s:%d", chroname, pos - pos % who.granularity)
	ky2 := fmt.Sprintf("%s:%d", chroname, pos - who.granularity - pos % who.granularity)

	_, OK := who.dict[ky1]
	if !OK{
		who.dict[ky1] = make(map[int] []string)
	}

	_, OK = who.dict[ky2]
	if !OK{
		who.dict[ky2] = make(map[int] []string)
	}

	_, OK = who.dict[ky1][pos]
	if !OK{
		who.dict[ky1][pos] = make([]string, 0)
	}

	_, OK = who.dict[ky2][pos]
	if !OK{
		who.dict[ky2][pos] = make([]string, 0)
	}


	who.dict[ky1][pos] = append(who.dict[ky1][pos] , value)
	who.dict[ky2][pos] = append(who.dict[ky2][pos] , value)
}


// this will lookup all events in the read, given that the reat starts at read_start_pos, and is shorter than granularity.
func (who * BucketedTable) Lookup(chroname string, read_start_pos, read_len int) (ret []string, found bool){
	ky := fmt.Sprintf("%s:%d", chroname, read_start_pos - read_start_pos % who.granularity)

	found = false
	poses, OK := who.dict[ky]
	if !OK{
		return
	}

	ret = make([]string, 0)
	for event_pos, values := range poses{
		if event_pos >= read_start_pos && event_pos < read_start_pos + read_len {
			ret = append(ret, values...)
			found = true
		}
	}

	return
}
