package util

import ("os";"bufio";"strings";"bytes";"compress/gzip";"fmt";"io";"os/exec";"errors")
import "encoding/binary"

type Qfile struct{
	filePointer	* os.File
	fileReader	* bufio.Reader
	rawReadCloser	* io.ReadCloser
	lastline	string
	reuseLastLine	bool
	fileInfo	string
	chro_table	map[string]*ChromosomeInfo
}

func Qstream(raw_reader *io.ReadCloser) (qfp * Qfile, err error){
	qfp = new(Qfile)
	qfp.fileReader = bufio.NewReaderSize(* raw_reader, 8192)
	qfp.rawReadCloser = raw_reader
	qfp.reuseLastLine = false
	return
}

func ReadUInt64LittleEndianFp(fp *os.File) (ret uint64, err error){
	bbuf := make([]byte, 8)
	nc, err := fp.Read(bbuf)
	if err != nil || nc != 8{
		if err != nil{return}
		if nc < 8{
			err = errors.New("Bad length")
			return
		}
	}
	ret = binary.LittleEndian.Uint64(bbuf)
	return
}

func ReadUInt16LittleEndianFp(fp *os.File) (ret uint16, err error){
	bbuf := make([]byte, 2)
	nc, err := fp.Read(bbuf)
	if err != nil || nc != 2{
		if err != nil{return}
		if nc < 2{
			err = errors.New("Bad length")
			return
		}
	}
	ret = binary.LittleEndian.Uint16(bbuf)
	return
}

func ReadUInt32LittleEndianFp(fp *os.File) (ret uint32, err error){
	bbuf := make([]byte, 4)
	nc, err := fp.Read(bbuf)
	if err != nil || nc != 4{
		if err != nil{return}
		if nc < 4{
			err = errors.New("Bad length")
			return
		}
	}
	ret = binary.LittleEndian.Uint32(bbuf)
	return
}

func QopenGz(gzfile string) (qfp * Qfile, err error){
	qfp = new(Qfile)
	if "STDIN://" == gzfile{
		qfp.filePointer = os.Stdin
	}else{
		qfp.filePointer, err = os.Open(gzfile)
	}
	if(err!=nil){
		return
	}

	gzFilePointer,err := gzip.NewReader(qfp.filePointer)
	if(err!=nil){
		return
	}

	qfp.fileReader = bufio.NewReaderSize(gzFilePointer, 8192)
	return

}

func Qopen(filename string) (qfp * Qfile, err error){

	qfp = new(Qfile)
	if "STDIN://" == filename{
		qfp.filePointer = os.Stdin
	}else{
		qfp.filePointer, err = os.Open(filename)
	}
	if(err!=nil){
		return
	}

	qfp.chro_table = make(map[string]*ChromosomeInfo)
	qfp.fileReader = bufio.NewReaderSize(qfp.filePointer, 8192)
	return
}

func (qf * Qfile) SetFileInfo(info string){
	qf.fileInfo = info
}

func (qf * Qfile) ReloadLine(){
	qf.reuseLastLine = true
}


func (qf * Qfile) Rewind(){
	if qf.filePointer == nil{
		panic("Cannot rewind a Qfile that was created on an io.reader.")
	}
	qf.filePointer.Seek(0, os.SEEK_SET)
	qf.fileReader = bufio.NewReaderSize(qf.filePointer, 8192)
}

func (qf * Qfile) TrimLine() (str string, err error){
	str, err = qf.Line()
	if err != nil{ return }
	str = strings.TrimSpace(str)
	return
}


func (qf * Qfile) Line() (str string, err error){
	if qf.reuseLastLine {
		qf.reuseLastLine = false
		str = qf.lastline
		return
	}

	str, err = qf.fileReader.ReadString('\n')
        if len(str)>0 {err=nil}
	qf.lastline = str
	return
}


func (qf * Qfile) Array() (strarr []string, err error){
	strarr, err = qf.Split("\t")
	return
}
func (qf * Qfile) Split(sep string) (strarr []string, err error){
	str, err := qf.Line()
	if(str!=""){
		strarr = strings.Split(strings.TrimSpace(str), sep)
		err = nil
	}
	return
}

func (qf * Qfile) Close() (err error){
	if qf.rawReadCloser != nil{
		(*qf.rawReadCloser).Close()
	}
	if os.Stdin != qf.filePointer && qf.filePointer != nil{
		err = qf.filePointer.Close()
	}
	return
}

func LoadTab(fn string, oldmap *map[string] string) map[string] string{
        ret := make(map[string] string, 10000)

	if oldmap != nil{
		for k,v := range (*oldmap){
			ret[k]=v
		}
	}
        fp, _:=Qopen (fn)
        for{
                farr, err := fp.Array()
		if err!=nil {
			break
		}
                ret[farr[0]]=farr[1]
        }

	fp.Close()
        return ret
}

func ReadFastA(fname string) (seqlist []string, seqmap map[string]string, err error){

	seqlist = make([]string,0)
	fp, err:=Qopen(fname)

	if err == nil{
		defer fp.Close()
		seqmap = make(map[string]string)
		var faseq bytes.Buffer
		seqname := ""
		for{
			fl, xerr := fp.Line()
			if xerr != nil || fl[0]=='>'{
				if seqname!=""{
					seqmap[seqname] = faseq.String()
					seqlist = append(seqlist, seqname)
				}

				faseq.Reset()

				if xerr==nil{
					seqname = strings.TrimSpace(fl[1:])
					if strings.Contains(seqname, " "){
						seqname = strings.Split(seqname, " ")[0]
					}
				}
			}else if xerr==nil{
				faseq.WriteString(strings.TrimSpace(fl))
			}

			if xerr!=nil{break}
		}
	}

	return
}

func ReadOnerowFastA(fname string) (rname, seq string, err error){
	fp, err:=Qopen(fname)

	seq = ""
	rname = ""
	rno := 0
	var faseq bytes.Buffer
	if err == nil{
		defer fp.Close()
		for{
			fl, xerr := fp.Line()
			if xerr!=nil{break}

			if rno == 0{
				rname = strings.TrimSpace(fl[1:])
			}else{
				faseq.WriteString(strings.TrimSpace(fl))
			}
			rno++
		}
		seq =faseq.String()
	}
	return
}

func IsFile(fname string) bool{
	if _, err := os.Stat(fname); err != nil {
		if os.IsNotExist(err) { return false }
	}
	return true
}


func WriteFastA(fp io.Writer, name , seq string){
	 WriteFastAEx(fp, name, seq, 70)
}

func WriteFastAEx(fp io.Writer, name , seq string, fasta_line_width int){
	wcur := 0
	fmt.Fprintf(fp, ">%s\n", name)


	for wcur < len(seq) {
		wcur_end := wcur + fasta_line_width
		wcur_end = MinInt(wcur_end, len(seq))
		str70 := seq[wcur:wcur_end]
		fmt.Fprintf(fp, "%s\n", str70)
		wcur += fasta_line_width
	}

}

func ListDir(path string)(files []string, err error){
	dirfp, err := os.Open(path)
	if err != nil{return}
	files, err = dirfp.Readdirnames(0)
	dirfp.Close()
	return
}

func TempFileName() string{
	return "/tmp/del4-YangLiao-"+ULID()
}

func DelTempFile(fn string){
	if strings.HasPrefix(fn,"/tmp/del4-YangLiao-"){
		os.Remove(fn)
	}
}

func Mv(o,n string){
	cmd := exec.Command("mv", "-f", o, n)
	cmd.Run()
}
