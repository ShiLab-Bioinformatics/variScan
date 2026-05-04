package util

import "encoding/binary"
import "os"
import "errors"
import "fmt"
import "strings"
import "path"

func load_bin_content(fx string) []int32{
        fobj,_ := os.Open(fx)
        finfo,_ := fobj.Stat()

        fsize := finfo.Size()
        if fsize < 4{
                panic("Unknown size!")
        }

        var bases int64
        bases = fsize/4

   //     fmt.Printf("Length=%d\n", bases)

	ret := make([]int32 , bases)
	for i := (int64)(0); i < bases ; {
		onesize := MinInt(10000, (int)(bases-i))
		var retone = make([]int32, onesize)
		err := binary.Read(fobj,  binary.LittleEndian, &retone)

		if err!=nil {
			panic(err)
		}

		for j := (int)(i); j< onesize + (int)(i); j++{
			ret[j] = retone[j-(int)(i)]
		}
		i += (int64)(onesize)
	}

        fobj.Close()

        return ret
}

func LoadCoverage(prefix string) (mc *MappingCoverage, err error){

	pref_files_ar := strings.Split(prefix, "/")
	pref_files := pref_files_ar[len(pref_files_ar)-1]
        pref_dir := path.Dir(prefix)
        dirfp, _ := os.Open(pref_dir)
        pref_names, _ := dirfp.Readdirnames(0)
        ret := make(map[string] []int32, 100)

        for _, nn := range pref_names{
		if strings.HasSuffix(nn, ".bin") && strings.HasPrefix(nn, pref_files){
			chro := strings.Split(nn,"--")[1]
			chro = strings.Split(chro, ".")[0]
		//	fmt.Printf("DIR=%s; chro:=%s\n" , pref_dir+"/"+nn, chro)
			ret[chro] = load_bin_content(pref_dir+"/"+nn)
		}
        }
        dirfp.Close()

	if len(ret)<1{
		err = errors.New("Unable to find any coverage bins for "+prefix)
	}else{
		mc = new(MappingCoverage)
		mc.binMap = ret
	}
	return
}

func (mc *MappingCoverage) Coverage(chro string, pos int) (coverage int, err error){
	coverage, err = mc.MeanCoverage(chro, pos, 1)
	return
}

func (mc *MappingCoverage) MeanCoverage(chro string, pos, bases int) ( meancoverage int, err error ){
	bins, found := mc.binMap[chro]
	if !found{
		err = errors.New("Unable to find chro: "+chro)
		return
	}
	if pos <0 || pos + bases > len(bins){
		err = errors.New(fmt.Sprintf("The range of request (%d ~ %d) excessed the range of chro: 0 ~ %d", pos, pos+bases, len(bins)) )
	}

	subv := 0
	for i := pos; i < pos+bases; i++{
		subv += int(bins[i])
	}

	meancoverage = subv / bases
	return
}

type MappingCoverage struct{
	binMap map[string] []int32
}
