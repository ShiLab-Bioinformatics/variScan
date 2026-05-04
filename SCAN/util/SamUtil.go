package util

import "strconv"
import "fmt"
import "errors"
import "unicode"
import "math/rand"
import "math"
import "bytes"
import "strings"

type ChromosomeInfo struct{
	Name string
	Length int
}

type SamBamReader interface{
	SAMRecordEx() (read_name string, flags int, chro string, pos , mapq int, cigar, mate_chro string, mate_pos, tlen int, seq, qual string, extra_columns []SAMAppendix, err error)
	SAMRecord() (read_name string, flags int, chro string, pos , mapq int, cigar, mate_chro string, mate_pos, tlen int, seq, qual string, err error)
	GetChromosomeInfo( chro_name string ) *ChromosomeInfo
}

type SAMAppendix struct{
	TagName string
	TagType	byte
	TagValue string
}

func (qf * Qfile) GetChromosomeInfo( chro_name string ) *ChromosomeInfo{
	return qf.chro_table[chro_name]
}

func (qf * Qfile) SAMRecord() (read_name string, flags int, chro string, pos , mapq int, cigar, mate_chro string, mate_pos, tlen int, seq, qual string, err error){
	read_name, flags, chro, pos, mapq, cigar, mate_chro, mate_pos, tlen, seq, qual, _, err = qf.SAMRecordEx()
	return
}

func ChroPos2Linear(chro string, pos int)(rpos int64, err error){
	err = nil
	if ! strings.HasPrefix(chro, "chr"){
		pos = 0
		err = errors.New("Not an ordinary chromosome")
		return
	}
	bpos := 0
	if unicode.IsDigit(rune(chro[3])){
		bpos, _ = strconv.Atoi(chro[3:])
	}else if chro[3] == 'X'{
		bpos = 31
	}else if chro[3] == 'Y'{
		bpos = 32
	}else if chro[3] == 'M'{
		bpos = 33
	}else{
		err = errors.New("Not an ordinary chromosome suffix")
	}
	rpos = int64(bpos) * int64(1000000000) + int64(pos)
	return
}

func SAMRecordParser(fli []string) (read_name string, flags int, chro string, pos , mapq int, cigar, mate_chro string, mate_pos, tlen int, seq, qual string, extra_columns []SAMAppendix, err error){
	read_name = fli[0]
	flags, _ = strconv.Atoi(fli[1])
	chro = fli[2]
	pos, _ = strconv.Atoi(fli[3])
	mapq, _ = strconv.Atoi(fli[4])
	cigar = fli[5]
	mate_chro = fli[6]
	mate_pos, _ = strconv.Atoi(fli[7])
	tlen , _ = strconv.Atoi(fli[8])
	seq = fli[9]
	qual = fli[10]

	extra_columns = make([]SAMAppendix, 0)
	for _,extra := range fli[11:]{
		exi := strings.Split(extra,":")
		if len(exi) == 3 && len(exi[0]) == 2 && len(exi[2])>0{
			new_appendix := SAMAppendix{exi[0], exi[1][0], exi[2]}
			extra_columns = append(extra_columns, new_appendix)
		}
	}
	return
}
func (qf * Qfile) SAMRecordEx() (read_name string, flags int, chro string, pos , mapq int, cigar, mate_chro string, mate_pos, tlen int, seq, qual string, extra_columns []SAMAppendix, err error){
	var fli []string
	for{
		fli, err = qf.Array()
		if err!= nil{
			return
		}
		if fli[0][0]!='@'{
			break
		}

		if fli[0]=="@SQ"{
			chro_name := ""
			chro_len := -1
			for _,tagstr := range fli[1:]{
				if strings.HasPrefix(tagstr,"SN:"){
					chro_name = tagstr[3:]
				} else if strings.HasPrefix(tagstr,"LN:"){
					chro_len , _ = strconv.Atoi(tagstr[3:])
				}
			}
			if chro_len>0 { qf.chro_table[chro_name] = & ChromosomeInfo{chro_name, chro_len} }
		}
	}
	//println("LLL0=",JoinStr("\t",fli))
	read_name, flags, chro, pos, mapq, cigar, mate_chro, mate_pos ,tlen, seq, qual, extra_columns, err = SAMRecordParser(fli)
	return
}

func ExIntTag(extra_columns []SAMAppendix , key string) int{
	for _,ex := range extra_columns{
		if ex.TagName==key{
			ret, ee := strconv.Atoi(ex.TagValue)
			if ee == nil{
				return ret
			}else{
				return -1;
			}
		}
	}
	return -1;
}

func GetIntTag(fl, tag string) int{
	ret := -1

	tag_str := fmt.Sprintf("\t%s:i:", tag)

	if strings.Contains(fl, tag_str){
		sli := strings.Split(fl, tag_str)
		sli = strings.Split(sli[1], "\t")
		var err error
		ret, err = strconv.Atoi(strings.TrimSpace(sli[0]))
		if err != nil{
			panic(err)
		}
	}

	return ret
}


func reverseBase(c byte) string{
	switch c{
		case 'A':
			return "T"
		case 'G':
			return "C"
		case 'C':
			return "G"
		case 'T':
			return "A"
	}
	return "N"
}

func SequenceErrorByRandomQualityStringEx(read string, error_scalling float64, refqual string) (newread, qual string, seq_error[]bool){
	var ATGC = [...]byte{'A','T','G','C'}
	var qprob []float64
	qual, qprob = SelectQualString(len(read), refqual)
	var buf bytes.Buffer
	seq_error = make([]bool, 0)

	for qpi,qp := range qprob{
		err := false
		nch := read[qpi]
		if rand.Float64() < qp * error_scalling * 1.3333333{
			och := nch
			nch = ATGC[int(rand.Float64() * 4)]
			err = och != nch
		}
		buf.WriteByte(nch)
		seq_error = append(seq_error, err)
	}
	newread = buf.String()
	return
}
func SequenceErrorByRandomQualityString(read string, error_scalling float64) (newread, qual string, seq_error[]bool){
	newread, qual, seq_error = SequenceErrorByRandomQualityStringEx(read, error_scalling, "SEQC-A")
	return
}

func SequenceErrorByQualityString(read, qstr string, error_scalling float64) (newread string, seq_error[]bool){
	var ATGC = [...]byte{'A','T','G','C'}
	var buf bytes.Buffer
	for qpi,qphred := range qstr{
		qp := math.Pow(10.0, -float64(qphred - quality_string_file_phred)/10.0)
		err := false
		nch := read[qpi]
		if rand.Float64() < qp * error_scalling * 1.3333333{
			och := nch
			nch = ATGC[int(rand.Float64() * 4)]
			err = och != nch
			//fmt.Printf("change: %c -> %c : %t\n", och, nch, err)
		}
		buf.WriteByte(nch)
		seq_error = append(seq_error, err)
	}
	newread = buf.String()
	return
}

func ReverseRead(bases, qual string) (rbases, rqual string){

	if len(bases) > 1000{
		var faseq bytes.Buffer
		for i:=len(bases) - 1 ; i >= 0; i--{
			faseq.WriteString(reverseBase(bases[i]))
		}

		rbases = faseq.String()

		if len(qual)>0{
			faseq.Reset()
			for i:=len(qual) - 1 ; i >= 0; i--{
				faseq.WriteString(string(qual[i]))
			}
			rqual = faseq.String()
		}
	}else{
		rbases = ""
		if len(bases)>0{
			for i:=0; i<len(bases); i++{
				rbases = reverseBase(bases[i]) + rbases
			}
		}

		if len(qual) > 0{
			rqual = ""
			if len(qual)>0{
				for i:=0; i<len(qual); i++{
					rqual = qual[i:i+1] + rqual
				}
			}
		}
	}

	return
}

const quality_string_file_phred = 33
var qualList map[string][]string

func load_qual_strs(reflib string)[]string{
	ret := make([]string, 0)
	filename := ""

	if reflib == "scRNA-136-R1"{
		filename="/home/liao_y/prj/Common/Libs/scripts/share/Nicolas-136-R1-qual.txt"
	}
	if reflib == "scRNA-136-R2"{
		filename="/home/liao_y/prj/Common/Libs/scripts/share/Nicolas-136-R2-qual.txt"
	}

	if reflib == "scRNA-Lisa-R1"{
		filename="/home/liao_y/prj/Common/Libs/scripts/share/4Org-Colon-R1-qual.txt"
	}
	if reflib == "scRNA-Lisa-R2"{
		filename="/home/liao_y/prj/Common/Libs/scripts/share/4Org-Colon-R2-qual.txt"
	}
	if filename == "" { panic("Unable to find a reflib for "+reflib) }

	fp, _ := Qopen(filename)

	for{
		ln, err:= fp.Array()
		if err!=nil{break}
		if len(ln)<1 || len(ln[0])<1{continue}
		ret=append(ret , ln[0])
	}

	fp.Close()
	return ret
}

func SelectQualString(rlen int, reflib string) ( rqual string, floats []float64){
	if nil == qualList{ qualList = make(map[string][]string) }
	qualstrs , OK := qualList[reflib]
	if ! OK{
		qualstrs = load_qual_strs(reflib)
		qualList[reflib] = qualstrs
	}

	qual := qualstrs[ rand.Intn(len(qualstrs))]

	sample_step := float64(len(qual))/float64(rlen) + 0.00000000001
	floats = make([]float64, 0)
	qual_new:=""
	for rli :=0 ; rli < rlen; rli++{
		qualidx := int(float64(rli) * sample_step)
		p_err := math.Pow(10.0, -float64(qual[qualidx] - quality_string_file_phred)/10.0)
		floats = append(floats,p_err)
		qual_new = qual_new + string(qual[qualidx] - quality_string_file_phred + 33)
	}
	rqual = qual_new
	return
}

