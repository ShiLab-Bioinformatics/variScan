package util

import "unicode"
import "fmt"
import "sort"
import "strconv"
import "strings"

type CigarOpt struct{
	OptType byte
	OptLen int
}

func ParseCigar(cigar string) []CigarOpt{
	if cigar == "*" { return nil }
	ret := make([]CigarOpt, 0)
	impt := 0
	for _, char  := range cigar{
		if unicode.IsDigit(char){
			impt = impt * 10 + int(char - '0')
		}else{
			var newpiece CigarOpt
			newpiece.OptType = byte(char)
			newpiece.OptLen = impt
			ret = append(ret,newpiece)

			impt = 0
		}
	}
	return ret
}

func Cigar14To13(cigar string) string{
	opts := ParseCigar(cigar)
	lastlen := 0
	ret := ""

	for _, opt := range opts{
		if opt.OptType=='=' || opt.OptType=='X' || opt.OptType=='M'{
			lastlen += opt.OptLen
		}else{
			if lastlen > 0{
				ret += fmt.Sprintf("%dM", lastlen)
				lastlen = 0
			}
			ret += fmt.Sprintf("%d%c", opt.OptLen, opt.OptType)
		}
	}
	if lastlen >0{
		ret += fmt.Sprintf("%dM", lastlen)
	}
	return ret
}


func ParseCigar14(cigar string) []CigarOpt{
	ret := make([]CigarOpt, 0)
	impt := 0
	last_len := 0
	last_opt := ' '
	for _, char := range cigar{
		if char == 'X'|| char=='=' { char = 'M' }

		if unicode.IsDigit(char){
			impt = impt * 10 + int(char - '0')
		}else{
			if char == last_opt{
				last_len += impt
			} else {
				if last_len > 0{
					var newpiece CigarOpt
					newpiece.OptType = byte(last_opt)
					newpiece.OptLen = last_len
					ret = append(ret,newpiece)
				}
				last_len = impt
				last_opt = char
			}

			impt = 0
		}
	}

	if last_len > 0{
		var newpiece CigarOpt
		newpiece.OptType = byte(last_opt)
		newpiece.OptLen = last_len
		ret = append(ret,newpiece)
	}

	if(false){
		for _,co := range ret{
			fmt.Printf("  Cigar: %c : %d\n", co.OptType, co.OptLen)
		}
	}
	return ret
}

func CigarIndelLen(opts []CigarOpt) int{
	ret := 0
	for _, co := range opts{
		if co.OptType == 'I' || co.OptType == 'D'{
			ret += co.OptLen
		}
	}

	return ret
}


func CigarChroLen(opts []CigarOpt) int{
	ret := 0
	for _, co := range opts{
		if co.OptType == 'M' || co.OptType == 'X' || co.OptType == '=' || co.OptType == 'D' || co.OptType == 'N' {
			ret += co.OptLen
		}
	}

	return ret
}



func CigarReadLen(opts []CigarOpt) int{
	ret := 0
	for _, co := range opts{
		if co.OptType == 'M' || co.OptType == 'X' || co.OptType == '=' || co.OptType == 'S' || co.OptType == 'I' || co.OptType == 'H'{
			ret += co.OptLen
		}
	}

	return ret
}


func CigarMappedLen(opts []CigarOpt, incI bool) int{
	ret := 0
	for _, co := range opts{
		if co.OptType == 'M' || co.OptType == 'X' || co.OptType == '=' || ( incI && co.OptType == 'I' ){
			ret += co.OptLen
		}
	}

	return ret
}

func FixSoftClipping(pos int, cigar string)(newpos int, newcigar string){
	opts := ParseCigar(cigar)
	first_S := 0
	last_S := 0
	last_M_index := 0

	for oii, opt := range opts{
		if opt.OptType == 'S'{
			if oii == 0{
				first_S = opt.OptLen
			}else{
				last_S = opt.OptLen
			}
		} else if opt.OptType == 'M'{
			last_M_index = oii
		}
	}

	is_first_M := true
	newcigar = ""
	for oii, opt := range opts{
		if opt.OptType == 'S'{
			continue
		} else if opt.OptType == 'M'{
			mlen := opt.OptLen

			if is_first_M{
				mlen += first_S
			}
			if oii == last_M_index{
				mlen += last_S
			}

			newcigar += fmt.Sprintf("%dM", mlen)
			is_first_M = false
		} else {
			newcigar += fmt.Sprintf("%d%c", opt.OptLen, opt.OptType)
		}
	}

	newpos = pos - first_S

	return
}

type MappedSection struct{
	Chro string
	Pos  int
	Cigar string
	ReadLength int
	IsNegative bool
	IsMainAlignment bool
	IsCoordinateGoingUp bool

	ConnectToLeft bool
	ConnectToRight bool

	ReadPosition int
	ChromosomeLength int
}

func CalcOverlap(fe GeneFeature, sec MappedSection) int{
	// GeneFeatures had the last base inclusive
	// MappedSection had the last base (Pos + ChromosomeLength) exclusive.
	if fe.Start + 1 <= sec.Pos + sec.ChromosomeLength{
		if fe.Stop >= sec.Pos{
			return MinInt(fe.Stop + 1, sec.Pos + sec.ChromosomeLength) - MaxInt(fe.Start, sec.Pos)
		}
	}
	return 0
}

func ParseCigarSections(chro string, pos int, cigar string, flags int) []MappedSection{
	ret := make([]MappedSection, 0)
	this_section_chro_start := pos
	this_cigar := ""
	this_read_length := 0
	this_section_read_start := 0
	this_chro_length := 0
	is_negative := (flags & 0x10)!=0

	cigar_opts  := ParseCigar(cigar)
	for copti, copt := range cigar_opts {
		if 'S' == copt.OptType{
			if copti < 1{
				this_section_read_start = copt.OptLen
			}
		}

		if 'M' == copt.OptType && copti == len(cigar_opts)-1{
			this_cigar += fmt.Sprintf("%dM", copt.OptLen)
			this_chro_length += copt.OptLen
			this_read_length += copt.OptLen
		}

		if 'N' == copt.OptType || 'D' == copt.OptType || copti == len(cigar_opts)-1 || 'I' == copt.OptType {
			nc := MappedSection{chro, this_section_chro_start, this_cigar, this_read_length, is_negative, true, true, true, true, this_section_read_start, this_chro_length}
			ret = append(ret, nc)

			this_section_read_start += this_read_length

			if 'I' == copt.OptType {
				this_section_chro_start += this_chro_length
				this_section_read_start += copt.OptLen
			} else if copti < len(cigar_opts)-1{
				this_section_chro_start += this_chro_length + copt.OptLen
			}

			this_read_length = 0
			this_chro_length = 0
			this_cigar = ""
		}

		if 'N' != copt.OptType  && 'S' != copt.OptType && 'D' != copt.OptType{
			this_cigar += fmt.Sprintf("%d%c", copt.OptLen, copt.OptType)
			if 'M' == copt.OptType{
				this_chro_length += copt.OptLen
				this_read_length += copt.OptLen
			}
		}
	}
	return ret
}

func ParseFusionSections(chro string, pos int, cigar string, flags int, appendix []SAMAppendix) []MappedSection{
	if flags & 4 == 4{ return make([]MappedSection, 0)}

	section_tab := make(map[int] MappedSection)
	cigar_opts  := ParseCigar(cigar)

	start_read_cursor := 0
	read_cursor := 0
	chromosomal_cursor := pos
	start_chromosomal_cursor := pos
	is_main_negative_strand := flags & 0x10 > 0
	is_second_read := flags & 0x80 > 0

	section_cigar := ""

	last_M_index := len(cigar_opts) -1
	for copti, copt := range cigar_opts {
		if copt.OptType == 'M'{
			last_M_index = copti
		}
	}

	for copti, copt := range cigar_opts {
		if copt.OptType != 'N' && copt.OptType != 'S'{
			section_cigar += fmt.Sprintf("%d%c", copt.OptLen, copt.OptType)
		}

		if copt.OptType == 'M'{
			read_cursor += copt.OptLen
			chromosomal_cursor += copt.OptLen
		}else if copt.OptType == 'S' || copt.OptType == 'I'{
			read_cursor += copt.OptLen
			if copt.OptType == 'S' {
				start_read_cursor = read_cursor
			}

		}else if copt.OptType == 'D'{
			chromosomal_cursor += copt.OptLen
		}

		if copt.OptType == 'N' || last_M_index == copti{
			/*if cigar=="5S115M"{
				fmt.Printf("RRR %d - %d -- %s at %s:%d\n", read_cursor, start_read_cursor, cigar, chro, pos)
			}*/
			section_tab[start_read_cursor] = MappedSection{chro, start_chromosomal_cursor, section_cigar, read_cursor - start_read_cursor, is_main_negative_strand != is_second_read, true, true, false, false, start_read_cursor, chromosomal_cursor - start_chromosomal_cursor}
			start_read_cursor = read_cursor
			section_cigar = ""

			if  copt.OptType == 'N'{
				chromosomal_cursor += copt.OptLen
				start_chromosomal_cursor = chromosomal_cursor
			}
		}

	}


	app_cigar := ""
	app_chro := ""
	app_pos := 0
	app_strand := ""
	if len(appendix) > 0{
		for _, appx := range appendix{
			if appx.TagName == "CG"{
				app_cigar = appx.TagValue
			}else if appx.TagName == "CP"{
				app_pos, _ = strconv.Atoi(appx.TagValue)
			}else if appx.TagName == "CT"{
				app_strand = appx.TagValue
			}else if appx.TagName == "CC"{
				app_chro = appx.TagValue
				//TODO: add into sections.

				is_apdx_negative_strand := "-" == app_strand
				apdx_cigar_opts := ParseCigar(app_cigar)
				read_cursor := 0
				if 'S' == apdx_cigar_opts[0].OptType{
					read_cursor = apdx_cigar_opts[0].OptLen
				}

				apdx_mapped_cigar := ""
				mapped_chro_len := 0
				mapped_read_len := 0
				for _, opt := range apdx_cigar_opts{
					if opt.OptType == 'M' || opt.OptType == 'I'{
						mapped_read_len += opt.OptLen
					}

					if opt.OptType == 'M' || opt.OptType == 'D'{
						mapped_chro_len += opt.OptLen
					}

					if opt.OptType != 'S'{
						apdx_mapped_cigar += fmt.Sprintf("%d%c", opt.OptLen, opt.OptType)
					}
				}

				if is_apdx_negative_strand == (is_main_negative_strand != is_second_read){
					section_tab[read_cursor] = MappedSection{app_chro, app_pos, apdx_mapped_cigar, mapped_read_len, is_apdx_negative_strand, false, true, false, false, read_cursor, mapped_chro_len}
				}else{
					section_tab[read_cursor] = MappedSection{app_chro, app_pos + mapped_chro_len - 1, apdx_mapped_cigar, mapped_read_len, is_apdx_negative_strand, false, false, false, false, read_cursor, mapped_chro_len}
				}
			}
		}
	}

	prev_encounter := 0
	next_encounter := 0
	double_false := false

	for rpos, pobj := range section_tab{

		has_prev := false
		for tpos, tobj := range section_tab{
			if tpos + tobj.ReadLength == rpos{ has_prev = true }
		}
		pobj.ConnectToLeft = has_prev

		_, has_next := section_tab[rpos + pobj.ReadLength]
		pobj.ConnectToRight = has_next

		if !has_prev{prev_encounter ++}
		if !has_next{next_encounter ++}
		if !has_prev && !has_next{double_false = true}

		section_tab[rpos] = pobj
	}

	if 1!= prev_encounter || 1 != next_encounter || (double_false && len(section_tab) > 1) || len(section_tab) <1{
		panic(fmt.Sprintf("Sections were not connected! Cigar=%s; Sections=%d", cigar, len(section_tab)))
	}

	sortints := make([]int,0)
	for rpos := range section_tab{
		sortints=append(sortints, rpos)
	}

	sort.Ints(sortints)

	ret := make([]MappedSection,0)

	for _, rpos:= range sortints{
		ret = append(ret, section_tab[rpos])
	}

	return ret
}


func LoadReportedJunctions(JuncOut_file string)( junctions []FusionEvent, err error){
	junction_fp, junction_err := Qopen(JuncOut_file)
	if junction_err != nil{ err = junction_err }else{
		junctions = make([]FusionEvent, 0)
		for{
			arrv, aerr := junction_fp.Array()
			if aerr != nil{ break }
			if arrv[0][0]=='#' {continue}

			ch1 := arrv[0]
			po1, _ := strconv.Atoi(arrv[1])
			ch2 := arrv[0]
			po2, _ := strconv.Atoi(arrv[2])
			nSup, _ := strconv.Atoi(arrv[4])

			pos_offsets:= strings.Split(arrv[10],",")

			pos_off1, _ := strconv.Atoi(pos_offsets[0])
			pos_off2, _ := strconv.Atoi(pos_offsets[1])

			is_cross := false
			ext_to_small_1 := true
			ext_to_small_2 := false

			ret := FusionEvent{ ch1, ch2, po1 + pos_off1, po2 - pos_off2 + 1, is_cross, ext_to_small_1, ext_to_small_2, 0, 0, "JUNCTION", nSup }

			junctions=append(junctions, ret)
		}

		junction_fp.Close()
	}
	return

}

func LoadFusions(ff string)( fusions []FusionEvent, err error){
	fusion_fp, fusion_err := Qopen(ff)
	if fusion_err != nil{ err = fusion_err }else{
		fusions = make([]FusionEvent, 0)
		for{
			arrv, aerr := fusion_fp.Array()
			if aerr != nil{ break }
			if arrv[0][0]=='#' {continue}
			if arrv[0]=="UNPAIRED"{
				ch1 := arrv[1]
				po1, _ := strconv.Atoi(arrv[2])
				ch2 := arrv[3]
				po2, _ := strconv.Atoi(arrv[4])
				nSup, _ := strconv.Atoi(arrv[6])

				is_cross := arrv[5]=="No"
				ext_to_small_1 := arrv[7]=="No"
				ext_to_small_2 := arrv[8]=="No"

				if(ext_to_small_1 == ext_to_small_2 && !is_cross){panic("Wrong strand cross 1")}
				if(ext_to_small_1 != ext_to_small_2 && is_cross){panic("Wrong strand cross 2")}

				ret := FusionEvent{ ch1, ch2, po1, po2, is_cross, ext_to_small_1, ext_to_small_2, 0, 0 , "FUSION", nSup}

				fusions=append(fusions, ret)
			}
		}

		fusion_fp.Close()
	}
	return
}

type FusionEvent struct{
	Chro1, Chro2  string
	Pos1,  Pos2   int
	Cross, ExtendToSmall1,  ExtendToSmall2  bool
	ExtensionLength1, ExtensionLength2 int
	EventType string
	NSup int
}

func (self FusionEvent) Key() string{
	k1 := fmt.Sprintf("%s:%d", self.Chro1, self.Pos1)
	k2 := fmt.Sprintf("%s:%d", self.Chro2, self.Pos2)
	if k1>k2{ return k2+"~"+k1}else{ return k1+"~"+k2}
}


func ExtractFusionEvents(sections []MappedSection, read_on_negative bool) []FusionEvent{
	ret := make([]FusionEvent, 0)
	if len(sections)<2{
		return ret
	}


	for ii, _ := range sections{
		if 0 == ii{continue}

		left := sections[ii-1]
		right := sections[ii]

		var extend_to_small_1 bool
		var p1 int

		if left.IsNegative == read_on_negative{
			p1 = left.Pos + left.ChromosomeLength - 1
			extend_to_small_1 = true
		}else{
			p1 = left.Pos - left.ChromosomeLength + 1
			extend_to_small_1 = false
		}

		extend_to_small_2 := right.IsNegative != read_on_negative
		p2 := right.Pos
		cross_connect := extend_to_small_2 == extend_to_small_1

		new_event := FusionEvent{left.Chro, right.Chro, p1, p2, cross_connect, extend_to_small_1, extend_to_small_2, left.ChromosomeLength, right.ChromosomeLength, "CIGAR", 1}
		ret = append(ret, new_event)
	}

	return ret
}

func ExtraTagInt(extra []SAMAppendix, tagname string) (found bool , tagvalue int){
	for _, ex := range extra{
		if ex.TagName == tagname{
			tagvalue, _ = strconv.Atoi(ex.TagValue)
			found = true
			return
		}
	}
	tagvalue = -1
	found = false
	return
}


// A should contain B; A should be much longer than B
func IsWellContaining(As, Ae, Bs, Be, min_hanging_len int) bool{

	if Bs > As + min_hanging_len && Be + min_hanging_len < Ae{
		return true
	}

	return false
}

// A and B should have at least min_hanging_len bases overlapping
func IsOverlapping(As, Ae, Bs, Be, min_hanging_len int) bool{
	max_start := MaxInt(As, Bs)
	min_end := MinInt(Ae, Be)

	return min_end - max_start >= min_hanging_len
}

func IsFusionRead(sections []MappedSection, max_intron int) bool{
	min_coor := 999999999
	max_coor := 0

	if len(sections)<1 {return false}
	chro1 := sections[0].Chro
	for _,sec := range sections{
		if chro1 != sec.Chro{ return true }
		min_coor = MinInt(min_coor, sec.Pos)
		max_coor = MaxInt(max_coor, sec.Pos + sec.ChromosomeLength)
	}

	return max_coor - min_coor > max_intron
}

func CompressCigar(cigar string)string{
	ret := ""
	optlen := -1
	sumlen := 0
	old_opt := 'x'
	for _, nch:= range cigar{
		if nch >= '0' && nch <= '9'{
			if optlen <0 { optlen = 0 }
			optlen = 10*optlen + int(nch - '0')
		}else{
			if optlen <0 { optlen = 1}
			if old_opt != nch && old_opt != 'x'{
				ret += fmt.Sprintf("%d%c", sumlen, old_opt)
				sumlen = 0
			}
			old_opt = nch
			sumlen += optlen
			optlen = -1
		}
	}
	if old_opt != 'x' && sumlen >0 { ret += fmt.Sprintf("%d%c", sumlen, old_opt) }
	return ret
}

func CropCigarOnChro(cigar string, first_base_loc, wanted_first_chro_base, wanted_last_chro_base int, add_S_before_after bool)(ret string, new_mapped_loc_after_S, skipped_rbases int){ // first_base_loc and wanted base have to be in the same system (either 0-based or 1-based)
	if wanted_last_chro_base <0 {wanted_last_chro_base = 0x7fffffff}

	ops := ParseCigar(cigar)
	chro_cursor := first_base_loc
	skipped_rbases =0
	last_kept_rbases :=-1
	firstS:=0
	ret = ""

	new_mapped_loc_after_S = -1
	read_cursor :=0
	for opi,op := range ops{
		//if op.OptType=='S' || op.OptType=='H'{ panic("No S and H should be included when we crop on chromosome coordinates.") }
		if op.OptType=='S' || op.OptType=='H'{
			if opi==0{ firstS = op.OptLen }
 			read_cursor += op.OptLen
		}
		if op.OptType=='I'{
			should_include_I := wanted_first_chro_base < chro_cursor && wanted_last_chro_base >= chro_cursor
			if should_include_I {
				ret+=fmt.Sprintf("%d%c", op.OptLen, op.OptType)
			}
			if wanted_first_chro_base >= chro_cursor {
				skipped_rbases += op.OptLen
			}
			read_cursor += op.OptLen
			if last_kept_rbases < read_cursor { last_kept_rbases = read_cursor -1}
		}
		if op.OptType=='M' || op.OptType=='D' || op.OptType=='N'{
			seg_output_start := -1
			seg_output_end := -1
			if wanted_first_chro_base <=chro_cursor + op.OptLen -1{
				seg_output_start = wanted_first_chro_base - chro_cursor
				if seg_output_start < 0{seg_output_start=0}
			}
			if wanted_last_chro_base >=chro_cursor{
				seg_output_end = wanted_last_chro_base - chro_cursor
				if seg_output_end > op.OptLen -1{seg_output_end=op.OptLen -1}
			}
			if seg_output_end >= 0 && seg_output_start >= 0{
				seg_used_len := seg_output_end - seg_output_start +1
				ret += fmt.Sprintf("%d%c", seg_used_len, op.OptType)
				if new_mapped_loc_after_S <0 {new_mapped_loc_after_S = chro_cursor + seg_output_start}
				if op.OptType== 'M' && last_kept_rbases < read_cursor + seg_output_end { last_kept_rbases = read_cursor + seg_output_end }
			}
			if seg_output_start <0 && op.OptType=='M'{ skipped_rbases += op.OptLen }
			if seg_output_start >0 && op.OptType=='M'{ skipped_rbases += seg_output_start }

			chro_cursor += op.OptLen
			if op.OptType=='M' {read_cursor += op.OptLen}
		}
	}
	if firstS >0 { ret = fmt.Sprintf("%dS", firstS)+ret }
	if add_S_before_after {
		if skipped_rbases>0{ ret = fmt.Sprintf("%dS", skipped_rbases)+ret}
		if last_kept_rbases>0 { ret = ret+fmt.Sprintf("%dS", read_cursor - (last_kept_rbases +1))}
	}
	ret = CompressCigar(ret)
	return
}


func CropCigarOnRead(cigar string, chro_pos_after_S, wanted_first_read_base_0based, wanted_last_read_base_0based int, add_S_before_after bool)(ret string, new_first_base_chro_loc, skipped_chrobases int){
	if wanted_last_read_base_0based <0 {wanted_last_read_base_0based = 0x7fffffff}

	ops := ParseCigar(cigar)
	read_cursor :=0
	skipped_chrobases = 0
	first_reported_rbase := 0x7fffffff
	last_kept_rbases :=-1
	ret = ""
	chro_cursor := chro_pos_after_S
	new_first_base_chro_loc =-1

	for _,op := range ops{
		if op.OptType=='D' || op.OptType=='N'{
			should_include_D := wanted_first_read_base_0based < read_cursor && wanted_last_read_base_0based >= read_cursor
			if should_include_D {
				ret+=fmt.Sprintf("%d%c", op.OptLen, op.OptType)
			}
			if wanted_first_read_base_0based >= read_cursor {
				skipped_chrobases += op.OptLen
                        }
			chro_cursor += op.OptLen
		}
		if op.OptType=='M' || op.OptType=='I' || op.OptType=='S' || op.OptType=='H'{
			seg_output_start := -1
			seg_output_end := -1
			if wanted_first_read_base_0based <=read_cursor + op.OptLen -1{
				seg_output_start = wanted_first_read_base_0based - read_cursor
				if seg_output_start <0 {seg_output_start=0}
			}
			if wanted_last_read_base_0based >=read_cursor{
				seg_output_end = wanted_last_read_base_0based - read_cursor
				if seg_output_end > op.OptLen -1{seg_output_end=op.OptLen -1}
			}
//if op.OptLen == 163{ println("INNERCHOP ", seg_output_start, seg_output_end, " CALC ",  wanted_first_read_base_0based, read_cursor, op.OptLen , "  >>",ret ) }
			if seg_output_end >= 0 && seg_output_start >= 0{
				seg_used_len := seg_output_end - seg_output_start +1
				ret += fmt.Sprintf("%d%c", seg_used_len, op.OptType)
				if last_kept_rbases < read_cursor + seg_output_end { last_kept_rbases = read_cursor + seg_output_end }
				if first_reported_rbase > read_cursor + seg_output_start { first_reported_rbase = read_cursor + seg_output_start }
				if op.OptType=='M' && new_first_base_chro_loc <0{ new_first_base_chro_loc = seg_output_start + chro_cursor }
			}

			if seg_output_start <0 && op.OptType=='M'{ skipped_chrobases += op.OptLen }
			if seg_output_start >0 && op.OptType=='M'{ skipped_chrobases += seg_output_start }
			if op.OptType=='M'{ chro_cursor += op.OptLen }
			read_cursor += op.OptLen
		}
	}
	if add_S_before_after {
		if first_reported_rbase>0 && first_reported_rbase<0x7fffffff{ ret = fmt.Sprintf("%dS", first_reported_rbase)+ret}
		if last_kept_rbases<read_cursor-1 { ret = ret+fmt.Sprintf("%dS", read_cursor - (last_kept_rbases+1))}
//println("INNER_BEFORE_COMPRESS", ret)
		ret = CompressCigar(ret)
	}
	return
}
