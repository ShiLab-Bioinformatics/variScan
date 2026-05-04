package util

import "os/exec"
import "bufio"
import "strings"
import "strconv"
import "io"
import "unicode"
import "encoding/binary"
import "os"
import "fmt"
import "bytes"
import "compress/gzip"
import "errors"

type RawBamChroInfo struct{
	ChroName string
	ChroLength int32
}

type ReadGroupInfo struct{
	GroupName string
}

type RawBAIChromosome struct{
	BinNumbers []int
	ChunksInBins [][][2]uint64
}

const BAI1_int = 0x01494142 // 'BAI\1' in Intel CPUs 

func GetBlockInfo(voff uint64) (bam_offset uint64, inblock_offset uint32 ){
	bam_offset = voff >>16
	inblock_offset = uint32(voff & 0xffff)
	return
}

func GetAllChromosomesBAI(fname string) (chroinfo []RawBAIChromosome, window16offsets [][]uint64, err error){
	var osfp * os.File
	osfp, err = os.Open(fname)
	if err != nil{return}

	magic_BAI1,_ := ReadUInt32LittleEndianFp(osfp);
	if(magic_BAI1 != BAI1_int){
		err = errors.New("No BAI1 magic int found! INT="+fmt.Sprintf("%d", magic_BAI1))
		return
	}

	chroinfo = make([]RawBAIChromosome,0)
	window16offsets = make([][]uint64,0)
	n_ref,_ := ReadUInt32LittleEndianFp(osfp)
	var ref_i, bin_i, chunk_i, win_i uint32
	for ref_i = 0; ref_i < n_ref; ref_i++{
		n_bin,_ := ReadUInt32LittleEndianFp(osfp)
		new_bai_chro := RawBAIChromosome{}
		new_bai_chro.BinNumbers = make([]int,0)
		new_bai_chro.ChunksInBins = make([][][2]uint64,0)

		for bin_i = 0; bin_i < n_bin; bin_i++{
			bin_number,_ := ReadUInt32LittleEndianFp(osfp)
			new_bai_chro.BinNumbers = append(new_bai_chro.BinNumbers, int(bin_number))
			n_chunk,_ := ReadUInt32LittleEndianFp(osfp)

			new_chunk_slice := make([][2]uint64,0)
			for chunk_i = 0; chunk_i < n_chunk; chunk_i++{
				chunk_start,_ := ReadUInt64LittleEndianFp(osfp)
				chunk_stop,_ := ReadUInt64LittleEndianFp(osfp)
				new_chunk_slice = append(new_chunk_slice, [2]uint64{chunk_start, chunk_stop})
			}
			new_bai_chro.ChunksInBins = append(new_bai_chro.ChunksInBins, new_chunk_slice)
		}
		chroinfo = append(chroinfo, new_bai_chro)
		n_windows,_ := ReadUInt32LittleEndianFp(osfp)
		new_wins := make([]uint64,0)
		for win_i =0; win_i < n_windows; win_i++{
			newi, _:= ReadUInt64LittleEndianFp(osfp)
			new_wins = append(new_wins, newi)
		}
		window16offsets = append(window16offsets, new_wins)
	}

	err = nil
	return

}

type RawBamFile struct{
	gzip_fp * gzip.Reader
	os_fp * os.File
	state int


	Chro_Info []RawBamChroInfo
	Read_Groups []ReadGroupInfo
}

func binReadHeader(fp * RawBamFile) error{
	if fp.state != 1{ panic("Shoudn't be called")}

	var txt_size int32
	err := fp.read_int_from_bam(&txt_size)
	if err != nil { return err }
	//println("HeaderText ",txt_size, "  ERR=", err)

	fp.Read_Groups = make([]ReadGroupInfo,0)
	headtxt_bin := make([]byte,txt_size)
	fp.Read(headtxt_bin)
	headtxt_buf := bytes.NewBuffer(headtxt_bin)
	headtxt_reader := bufio.NewReader(headtxt_buf)
	for{
		hline, herr :=headtxt_reader.ReadString('\n')
		if herr != nil{break}
		hline = strings.TrimSpace(hline)
		//fmt.Printf("%s\n", hline)
		if strings.HasPrefix(hline, "@RG\t"){
			secs := strings.Split(hline, "\t")
			for _, hsec := range secs[1:]{
				if strings.HasPrefix(hsec, "ID:"){
					grpname := hsec[3:]
					fp.Read_Groups = append(fp.Read_Groups, ReadGroupInfo{grpname})
					fmt.Printf("RG: `%s` (%d)\n", grpname, len(grpname))
				}
			}
		}
	}

	var ref_no int32
	err = fp.read_int_from_bam(&ref_no)
	if err != nil { return err }
	//println("RefSeqs ", ref_no, "  ERR=", err)

	fp.Chro_Info = make([]RawBamChroInfo,0)
	var ref_i int32
	for ref_i = 0; ref_i < ref_no; ref_i++{
		var refname_len, ref_seq_len int32
		err = fp.read_int_from_bam(&refname_len)
		if err != nil { return err }
		refname_bin := make([]byte , refname_len)
		rlen :=0
		rlen , err = fp.read_bytes_from_bam(refname_bin)
		if err != nil || int32(rlen) != refname_len{
			err = errors.New("Unable to load sequence name")
			return err
		}

		err = fp.read_int_from_bam(&ref_seq_len)
		if err != nil { return err }

		seqname := string(refname_bin[:len(refname_bin)-1])
		//println("REF_SEQ:",seqname, len(seqname), ref_seq_len)
		fp.Chro_Info=append(fp.Chro_Info, RawBamChroInfo{seqname, ref_seq_len})
	}
	fp.state = 10

	return nil
}

func RawBamOpen(fn string) ( fp * RawBamFile, err error ){
	fp = &RawBamFile{}
	os_FP, err := os.Open(fn)
	fp.state = -1
	fp.os_fp = os_FP 
	if err == nil{
		fp.gzip_fp,_ = gzip.NewReader(os_FP)
		bbuf := make([]byte, 4)
		rlen, Xerr := fp.read_bytes_from_bam(bbuf)
		if Xerr == nil && rlen == 4 {
			if bbuf[3]==1 && bbuf[0]=='B' && bbuf[1]=='A' && bbuf[2]=='M'{
				fp.state = 1
				Yerr := binReadHeader(fp)
				if Yerr != nil{
					err = Yerr
				}
			}else{
				err = errors.New("Unable to find 1BAM from header!")
			}
		}else{
			err = errors.New("Unable to find header!")
		}
	}
	return
}

func (rb *RawBamFile) Close(){
	if rb.state > 0{
		rb.os_fp.Close()
	}
}

func (rb *RawBamFile) read_int_from_bam(data interface{}) error{
	var ibinst [8]byte
	switch data := data.(type) {
		case *int32:
			var ibin []byte = ibinst [:4]
			Xrlen, Xerr := rb.read_bytes_from_bam(ibin)
			if Xrlen != 4{
				return errors.New("EOF.")
			}else if Xerr == nil{
				*data = int32(binary.LittleEndian.Uint32(ibin))
				return nil
			}else {
				return Xerr
			}
	}

	return errors.New("What a type??")
}

func (rb *RawBamFile) Read(data []byte) (rlen int, err error){
	rlen, err = rb.read_bytes_from_bam(data)
	return
}
func (rb *RawBamFile) read_bytes_from_bam(data []byte) (rlen int, err error){
	read_cur := 0
	for read_cur < len(data) {
		Xrlen, Xerr :=rb.gzip_fp.Read(data[read_cur:])
		XerrStr :=""
		if Xerr != nil{ XerrStr = Xerr.Error() }
		_ = XerrStr
		//println( "Read " ,Xrlen," at ",read_cur ," until " , cap(data), " Err =", XerrStr)
		if Xerr != nil{
			XXerr := rb.reopen_bam()
			if XXerr != nil{
				err = XXerr
				return
			}
		}
		read_cur += Xrlen
	}

	rlen = read_cur
	return
}

func (rb *RawBamFile) reopen_bam()error{
	rb.gzip_fp.Close()
	var err error
	rb.gzip_fp, err = gzip.NewReader(rb.os_fp)
	errstr :="NIL"
	if err != nil{ errstr = err.Error() }
	_=errstr
	//println("REOPEN_BAM: ", errstr)
	return err
}

func (rb *RawBamFile) NextReadBin()(ret []byte, vfile_pos_rstart int64, err error) {
	err = nil
	if rb.state <=0{
		err = errors.New("BAM fp invalid")
	}else if rb.state == 1{ // before reading header
		binReadHeader(rb)
	}else if rb.state == 10{ // reading alignments
		var binlen int32
		err = rb.read_int_from_bam(&binlen)
		if err != nil{
			err = errors.New("File terminated")
			return
		}
		ret = make([]byte, binlen)
		rlen :=0
		rlen, err = rb.read_bytes_from_bam(ret)
		if err != nil || int32(rlen) != binlen{
			err = errors.New("File terminated")
			return
		}
	}
	return
}

func (rb *RawBamFile) ReadBinBasicInfo(b []byte) ( readname string, flag int32, chro string, pos, cigar_opts, read_len int32, extra_data []byte){
	ref_id := int32(binary.LittleEndian.Uint32(b[0:4]))
	if ref_id < 0{ chro="*" } else { chro = rb.Chro_Info[ref_id].ChroName  }
	pos = int32(binary.LittleEndian.Uint32(b[4:8]))
	bin_mq_nl := int32(binary.LittleEndian.Uint32(b[8:12]))
	flag_nc  := int32(binary.LittleEndian.Uint32(b[12:16]))
	read_len = int32(binary.LittleEndian.Uint32(b[16:20]))
	flag = flag_nc >>16
	cigar_opts = flag_nc & 0xffff
	readname_len := bin_mq_nl&0xff
	readname = string( b[32:32+readname_len-1] )

	extra_pos := 32+readname_len+cigar_opts*4+read_len+(read_len+1)/2
	if extra_pos < int32(len(b)){
		extra_data = b[extra_pos:]
	}else{
		extra_data = make([]byte,0)
	}
	return
}

func next_tag_pos(b []byte, curs int) int{
	tpc := b[curs+2]
	datalen :=1
	if tpc == 'i' || tpc == 'I' || tpc=='f' { datalen = 4  }
	if tpc == 's' || tpc == 'S'{ datalen = 2  }
	if tpc == 'Z' || tpc == 'H'{
		string_start := curs+3
		for strcurs := string_start; strcurs < len(b); strcurs++{
			if b[strcurs]==0{
				datalen = strcurs - string_start +1
				break
			}
		}
	}

	//fmt.Printf("PASSTAG %c%c -- %c DLEN %d\n", b[curs], b[curs+1], tpc, datalen)
	return curs + 3 + datalen
}

func tag_name_type_match(b []byte, curs int,  tagname string, ctype rune) bool{
	if b[curs]==tagname[0] && b[curs+1]==tagname[1] && rune(b[curs+2])==ctype {return true}
	return false
}

func (rb *RawBamFile) FindAllStringTags(b []byte, tagname string) []string{
	curs := 0
	ret := make([]string, 0)
	for curs < len(b){
		//println("TEST_CURS: ", curs, "<",len(b))
		if tag_name_type_match(b, curs, tagname, 'Z'){
			string_start := curs+3
			for strcurs := string_start; strcurs < len(b); strcurs++{
				if b[strcurs]==0{
					ret = append(ret, string(b[string_start:strcurs]))
					break
				}
			}
		}
		curs = next_tag_pos(b, curs)
	}
	return ret
}

type BamReader struct{
	out_fp io.ReadCloser
	out_reader * bufio.Reader
	cmd * exec.Cmd
	lastLine, bam_header string
	chro_table map[string]*ChromosomeInfo
	reloaded bool
}

func BAMopen(fname string) (br *BamReader, err error){
	br = new(BamReader)
	//br.cmd = exec.Command("/usr/local/bioinfsoftware/samtools/samtools-0.1.18/bin/samtools", "view", fname)
	br.cmd = exec.Command("samtools", "view", "-h", fname)
	br.out_fp, err = br.cmd.StdoutPipe()
	br.chro_table = make(map[string]*ChromosomeInfo)
	if err != nil {
		return
	}

	br.out_reader = bufio.NewReader(br.out_fp)
	err = br.cmd.Start()
	if err == nil{
		for{
			hstr, herr := br.Line()
			if herr != nil || hstr[0]!='@'{
				br.lastLine = hstr
				br.reloaded = true
				break
			}
			br.bam_header += hstr
			if strings.HasPrefix(hstr, "@SQ\t"){
				hsi := strings.Split(strings.TrimSpace(hstr), "\t")
				chro_name := ""
				chro_len := -1
				for _,tagstr := range hsi[1:]{
					if strings.HasPrefix(tagstr,"SN:"){
						chro_name = tagstr[3:]
					} else if strings.HasPrefix(tagstr,"LN:"){
						//println("'"+tagstr[3:]+"'")
						chro_len , _ = strconv.Atoi(tagstr[3:])
					}
				}
				if chro_len>0 { br.chro_table[chro_name] = & ChromosomeInfo{chro_name, chro_len} }
			}
		}
	}
	return
}

func (br * BamReader) Close(){
	br.out_fp.Close()
	br.cmd.Wait()
	return
}

func (br * BamReader) Line() (st string, err error){
	if br.reloaded{
		st = br.lastLine
		br.reloaded = false
	} else {
		st, err = br.out_reader.ReadString('\n')
		if err == nil{ br.lastLine = st }
		if st == ""{ err = errors.New("EOF") }
	}
	return
}


func (br * BamReader) ReloadLine(){
	br.reloaded = true
}

func (qf * BamReader) GetChromosomeInfo(chro_name string) *ChromosomeInfo{
	return qf.chro_table[chro_name]
}

func (qf * BamReader) SAMRecord() (read_name string, flags int, chro string, pos , mapq int, cigar, mate_chro string, mate_pos, tlen int, seq, qual string, err error){
	read_name, flags, chro, pos, mapq, cigar, mate_chro, mate_pos, tlen, seq, qual, _, err = qf.SAMRecordEx()
	return
}

func (qf * BamReader) BamHeader()string{
	return qf.bam_header
}

func (qf * BamReader) SAMRecordEx() (read_name string, flags int, chro string, pos , mapq int, cigar, mate_chro string, mate_pos, tlen int, seq, qual string, extra_columns []SAMAppendix, err error){
	var fli []string
	for{
		fli, err = qf.Array()
		if err!= nil{
			return
		}
		if fli[0][0]!='@' && len( fli ) >9&& len(fli[0]) < 256 && len(fli[1]) < 6 && len(fli[3]) < 20 && len(fli[4]) < 4  && unicode.IsDigit(rune(fli[1][0])) && unicode.IsDigit(rune(fli[3][0]))  && unicode.IsDigit(rune(fli[4][0])){
			break
		}
	}
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
		if len(exi)<3{
			//fmt.Println( read_name , extra)
		}else{
			new_appendix := SAMAppendix{exi[0], exi[1][0], exi[2]}
			extra_columns = append(extra_columns, new_appendix)
		}
	}

	return
}


func (qf * BamReader) Array() (strarr []string, err error){
	str, err := qf.Line()
	if(err==nil){
		strarr = strings.Split(strings.TrimSpace(str), "\t")
	}
	return
}


