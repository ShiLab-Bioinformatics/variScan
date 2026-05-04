package main

import "github.com/ShiLab-Bioinformatics/variScan/SCAN/util"
import "flag"
import "crypto/sha1"
import "io"
import "strings"
import "encoding/binary"
import "strconv"
import "fmt"
import "bytes"
import "sort"
import "os"

var max_MM int
const verbose = false
const reportAllBest = true

func main(){
    var binf1, binf2, r1f, r2f, libf string
    var rlen, nchros int

    flag.StringVar(&libf,"lib","","Lib CSV file")
    flag.StringVar(&binf1,"binf1","","Binary file for R1")
    flag.StringVar(&binf2,"binf2","","Binary file for R2")

    flag.StringVar(&r1f,"R1","","fastq file for R1 (zcat it if necessary)")
    flag.StringVar(&r2f,"R2","","fastq file for R2 (zcat it if necessary)")

    flag.IntVar(&rlen,"rlen",0,"Maximum read length (must be the same as rlen for running find_best_align)")
    flag.IntVar(&max_MM,"maxMM",0,"Alloed max mismatch at any end (the other end can have many mismatches)")

    flag.Parse()

    libcsv,_ := util.Qopen(libf)
    seqnames := make(map[int]string)
    for {
        fli, err:= libcsv.Split(",")
        if err !=nil{break}

        seqnames[len(seqnames)+1]=fli[1]
    }
    libcsv.Close()
    nchros = len(seqnames)
    if verbose {println("NCHROS=",nchros)}

    unit_size := rlen + 8*nchros
    if unit_size<1{panic("No rlen nor chros are given")}

    var BinOffsets [2]map[[20]byte] int
    rbin := make([]byte, unit_size)

    binfps := make([]*os.File,2)
    for binno := 1; binno<=2; binno++{
        binf := binf1
        if binno==2{binf=binf2}
        binfps[binno-1],_ = os.Open(binf)

        offtab := make(map[[20]byte] int)
        fpos := 0
        for{
            qn,_ := io.ReadFull(binfps[binno-1], rbin)
            if qn < unit_size {break}
            rseq := string(rbin[0:rlen])
            hashv := sha1.Sum(rbin[0:rlen])
            offtab [hashv] = fpos
            fpos += qn
            if len(offtab) % 200000 == 0{ fmt.Printf("tabsize %d, seq=%s   hash=%x\n", len(offtab),rseq, hashv)}
            //if strings.Contains(rseq,"TTGGATTCTTTAATAACAGCTATAACTACAAATGGAGCTCATCCTAGTAAAATGGTTACCATACAGAGAACATTGGATGGGAGGCTTCAGGTGGCTGGTCGGAAAGGATTTCCTCATGTGATCTATGCCCGTCTCTGGAGGTGGCCTGAT"){ fmt.Printf("tabsize %d, seq=%s   hash=%x\n", len(offtab),rseq, hashv)}
     //       if fpos > 2*69666{break}
        }
        BinOffsets[binno-1] = offtab
    }

   fq1f, _ := util.Qopen(r1f)
   fq2f, _ := util.Qopen(r2f)
   nline :=0
   rname := ""
   for{
      if nline % 800000 ==0{ fmt.Printf("Run reads %d\n", nline/4) }

      q1seq,err := fq1f.Line()
      q2seq,_ := fq2f.Line()
      if err!= nil{break}
      q1seq = strings.TrimSpace(q1seq)
      q2seq = strings.TrimSpace(q2seq)
      if nline %4 ==0{
         rname  = strings.Split(q1seq[1:]," ")[0]
         rname2:= strings.Split(q2seq[1:]," ")[0]
         if rname!=rname2{panic("RNAMES MISMA R1R2")}
      }else if nline %4 == 1{
         high_total_mm := 9999999
         high_score := -1
         high_occurance := 0
         high_r1_refstart :=-99999
         high_r2_refstart :=-99999
         mm_R1:= 9999
         mm_R2:= 9999
         var high_chros []int

         r1_matches := make(map[int][3]int)
         r2_matches := make(map[int][3]int)
         libnos_sorted := make([]int,0)
         for ri := 1; ri <=2; ri++{
            rseq := q1seq
            if ri==2{rseq = q2seq}
            qbyte := []byte(rseq)
            if len(qbyte) < rlen{
               padding := bytes.Repeat([]byte{0x20}, rlen-len(qbyte))
               qbyte = append(qbyte, padding...)
            }
            hashv := sha1.Sum(qbyte)
            fpos, OK := BinOffsets[ri-1][hashv]
            if ! OK{panic("A read isn't found "+rseq +"   Rno  "+strconv.Itoa(ri))}

            binfps[ri-1].Seek(int64(fpos),0)
            qn,_ := io.ReadFull(binfps[ri-1], rbin)
            if qn != unit_size {panic("Unable to read a readbin.")}
            read_hash := sha1.Sum(rbin[0:rlen])
            if read_hash != hashv {panic("Wrong value in a readbin.")}

            for libnoi := 0; libnoi <nchros; libnoi++ {
               qlpos := rlen+libnoi*8
               libno := int(binary.LittleEndian.Uint16(rbin[qlpos:qlpos+2]))

               qlpos = rlen+libnoi*8+4
               bin_match := int(binary.LittleEndian.Uint16(rbin[qlpos:qlpos+2]))

               bin_start := 0
               bin_mm := 9999

               if bin_match>0{
                  if bin_match>rlen{panic("DATA FORMAT WRONG!")}
                  qlpos = rlen+libnoi*8+2
                  bin_start = int(binary.LittleEndian.Uint16(rbin[qlpos:qlpos+2]))
                  qlpos = rlen+libnoi*8+6
                  bin_mm = int(binary.LittleEndian.Uint16(rbin[qlpos:qlpos+2]))
               }
               rxm3 := [3]int{bin_match,bin_start,bin_mm}
               if ri == 2{ r2_matches[libno] = rxm3 } else {
                  r1_matches[libno] = rxm3
                  libnos_sorted = append(libnos_sorted,libno)
              }
            }
         }
         sort.Ints(libnos_sorted)
         for _,seq1no := range libnos_sorted{
            score1v := r1_matches[seq1no]
            score2v := r2_matches[seq1no]
            if score1v[2] > max_MM && score2v[2] > max_MM{ continue }
            if score1v[2]>1000 || score2v[2]>1000 {continue}
            if score1v[0]<1 || score2v[0]<1 {continue}
            score1 := score1v[0]
            score2 := score2v[0]
            my_score := score2 + score1
            if my_score > high_score || ( my_score == high_score && high_total_mm > score1v[2]+score2v[2]){
               high_chros = make([]int,1)
               high_chros[0] = seq1no
               high_score = my_score
               high_r1_refstart = score1v[1] // only works for the single-best mode.
               high_r2_refstart = score2v[1] // only works for the single-best mode.
               high_occurance = 1
               high_total_mm = score1v[2]+score2v[2]
               mm_R1 = score1v[2]
               mm_R2 = score2v[2]
            }else if my_score == high_score && high_total_mm == score1v[2]+score2v[2]{
               if reportAllBest {high_chros = append(high_chros, seq1no) }
               high_occurance ++
            }
         }

         chro1 := "NA"
         for outi := 0; outi <2; outi++{
            if len(high_chros)<=outi{ break }
            chro1 = seqnames[high_chros[outi]]
            if high_r1_refstart > 32767{ high_r1_refstart = high_r1_refstart-65536 }
            if high_r2_refstart > 32767{ high_r2_refstart = high_r2_refstart-65536 }
            fmt.Printf("READ %s has_score %d mapped_loc_R1_R2 %d %d rlen %d total_num_mapped_sequences %d mapped_to %s num_MM_R1_R2 %d %d [%s]\n", rname, high_score, high_r1_refstart, high_r2_refstart, rlen, high_occurance, strings.ReplaceAll(chro1," ","_"), mm_R1, mm_R2, "ONLY_POSITIVE_CONSIDERED")
         }
         if len(high_chros) <1{ fmt.Printf("READ %s has NO_results\n", rname) }
      }
      nline++
   }
   binfps[0].Close()
   binfps[1].Close()
   fmt.Println("ALL_MATCHING_TASK_FINISHED")
}
