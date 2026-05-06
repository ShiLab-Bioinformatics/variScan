package main

import "github.com/ShiLab-Bioinformatics/variScan/SCAN/util"
import "flag"
import "fmt"
import "strings"
import "os"
import "bufio"
import "time"
import "bytes"
import "encoding/binary"

const verbose=false
var workers int
var max_read_len int

type chanRes struct{
  rname string
  high_chros []int
  high_score, high_occurance, high_r1_refstart, high_r2_refstart,rlen, high_mm, mmR1, mmR2 int
  high_is_neg bool
}

// [Matched/Mismatched] -> [ACGT] -> ref_seq_base_no -> list_chro_no
type alike_seq_index struct{
   seqs map[int]string
   ints [2][4][][]int
   shrotest_refseq_len int
}

func make_alike_seq_index(csv *util.Qfile, removeRepeatSeq bool)(ret alike_seq_index, nseq int, refseq1 string, no_chrs[]bool, name_chros map[int]string){
  refseqs := make(map[int]string)
  seqreverse := make(map[string]int)
  name_chros = make(map[int]string)
  max_seqlen := 0
  min_seqlen := 999999999
  nseq = 0
  no_chrs = make([]bool, 1, 1)
  no_chrs [0] = true // chr_no starts 1
  ret.seqs = make(map[int]string,0)
  for{
    qline, err := csv.Line()
    if err != nil {break}
    if len(qline)<10{continue}
    qline = strings.TrimSpace(qline)

    nseq=1+nseq
    qnames := strings.Split(qline,",")
    name_chros[nseq]=qnames[1]

    if removeRepeatSeq{
      old_chro_No, OK := seqreverse[qnames[0]]
      if OK{
        delete(refseqs,old_chro_No)
        no_chrs[old_chro_No]=true
      }
    }
    refseqs[nseq] = qnames[0]
    seqreverse[qnames[0]] = nseq
    no_chrs = append(no_chrs, false)
    if max_seqlen < len(qnames[0]) { max_seqlen = len(qnames[0]) }
    if min_seqlen > len(qnames[0]) { min_seqlen = len(qnames[0]) }
  }
  if verbose{println("NUM_HEAD_SEQ",nseq)}
  ret.shrotest_refseq_len = min_seqlen

  for rii:=0; rii < 19999; rii++{
    rseq, OK := refseqs[rii]
    if OK && len(rseq)<max_seqlen{
        ret.seqs[rii] = rseq
    }
    if OK && rii > 11111{ panic("Not expected to have very many ref seqs!") }
  }


  for atgci :=0; atgci<4; atgci++{
    ret.ints[0][atgci] = make([][]int,max_seqlen,max_seqlen)
    ret.ints[1][atgci] = make([][]int,max_seqlen,max_seqlen)
  }

  for atgci :=0; atgci<4; atgci++{
    for base_i :=0 ; base_i < max_seqlen ; base_i++{
      ret.ints[0][atgci][base_i] = make([]int,0)
      ret.ints[1][atgci][base_i] = make([]int,0)
    }
  }

  for seq_i :=1; seq_i <=nseq; seq_i++{
    seq, OK := refseqs[seq_i]
    if ! OK || len(seq)!=max_seqlen{continue}
    for base_i :=0 ; base_i < len(seq) ; base_i++{
      basech := seq[base_i]
      atgc_i := 0
      if basech == 'C'{atgc_i=1}
      if basech == 'G'{atgc_i=2}
      if basech == 'T'{atgc_i=3}
      if basech == 'N'{atgc_i=-999}

      if atgc_i >= 0{
         ret.ints[0][atgc_i][base_i]=append(ret.ints[0][atgc_i][base_i], seq_i)
         for noX:=0; noX<=3; noX++{
            if noX!=atgc_i { ret.ints[1][noX][base_i]=append(ret.ints[1][noX][base_i], seq_i) }
         }
      }
    }
  }
  refseq1 = refseqs[1]
  return
}


func main(){
  var old_result, R1_fastq, lib_fasta, ofname string
  var removeRepeatSeq, is_R2 bool

  flag.StringVar(&ofname,"outfile","","Outfile (binary) name.")
  flag.StringVar(&old_result,"oldres","","Result file containing reads that don't need to process again.")
  flag.StringVar(&R1_fastq,"R1","","Read 1 FASTQ file name. Must be zcatted if it is a gz file.")
  flag.StringVar(&lib_fasta,"lib","","CSV file containing reference sequences (ref_seq,name).")
  flag.BoolVar(&removeRepeatSeq,"removerep",false,"Remove repeating sequences from the CSV library reference file. Only keep the last copy in file.")
  flag.BoolVar(&is_R2,"modeR2",false,"The input reads are from R2 fastq (i.e., needs reverse for comparison).")
  flag.IntVar(&workers,"threads",8,"Number of threads for running the alignment.")
  flag.IntVar(&max_read_len,"rlen",-1,"Maximum read length in FASTQ files.")
  flag.Parse()

  if max_read_len <1{panic("Maximum read length must be specified!")}
  old_res_map := make(map[string]bool)
  Qold, cnterr := util.Qopen(old_result)
  if cnterr == nil{
    for{
      fl, err := Qold.Line()
      if err != nil{break}
      fli := strings.Split(fl," ")
      if fli[0]=="READ"{
        hadname := strings.Clone(fli[1])
        old_res_map[hadname]=true
      }
    }
  }

  r1_fp, err1 := util.Qopen(R1_fastq)
  lib_fp, err3 := util.Qopen(lib_fasta)

  is_wrong := false

  if err1 != nil{
    println("Error: R1 FASTQ file not found.")
    is_wrong = true
  }
  if is_wrong == false && err3 != nil{
    println("Error: reference FASTA file not found.")
    is_wrong = true
  }

  ofp,_:= os.Create(ofname)
  writer := bufio.NewWriterSize(ofp, 1024*1024*10)
  if is_wrong == false{ run_job(r1_fp, lib_fp,removeRepeatSeq, is_R2, writer) }
  writer.Flush()
  ofp.Close()
}

func run_job(r1, lib *util.Qfile, removeRepeatSeq, is_R2 bool, ofp  *bufio.Writer){
  alike_index, max_lib_no, refseq1, no_chrs, _ := make_alike_seq_index(lib, removeRepeatSeq)

  novalid_chrs :=0
  for _, nov := range no_chrs{
     if nov {novalid_chrs ++}
  }

  no_seqs_for_seqs := make([]int,0)
  for nni := 0; nni <= 19999; nni++{
     _,OK := alike_index.seqs[nni]
     if OK{ no_seqs_for_seqs=append(no_seqs_for_seqs,nni) }
     if OK && nni>11111{ panic("Lib is too large.") }
  }

  inchs := make([]chan jobdef,workers)
  outchs := make([]chan jobres,workers)
  for wii := 0; wii < workers; wii++{
    inchs[wii] = make(chan jobdef)
    outchs[wii] = make(chan jobres, 1)
    go worker(inchs[wii], outchs[wii], refseq1, max_lib_no, alike_index, novalid_chrs, no_chrs, no_seqs_for_seqs)
  }

  my_worker:=0
  run_no := 0
  end_workers:= 0
  fl0:=make([]string,workers)
  for{
    if my_worker == workers{ my_worker=0 }

    fl , err:= r1.Line()
    newfl:=""
    tps := time.Now()
    if err == nil{
       newfl = strings.Split(strings.TrimSpace(fl) ," ")[1]
       if max_read_len < len(newfl){panic("Maximum read length is not sufficient!")}
       if alike_index.shrotest_refseq_len <= len(newfl){panic("A read has the same length or is longer than the shortest reference sequence!")}
       qseq := newfl
       if is_R2 { qseq,_ = util.ReverseRead(qseq,"") }
       inchs[my_worker] <- jobdef(qseq)
    }else{ end_workers++ }
    if end_workers == 1+workers{break}

    run_no ++
    if run_no <= workers{
       fl0[my_worker]=newfl
       my_worker++
       continue
    }

    jrr := <-outchs[my_worker]

    tp0 := time.Now()
    add_one_rec(ofp, strings.TrimSpace(fl0[my_worker]))
    fl0[my_worker]=newfl
    nwt := 0
    for libseqnum, n_matched := range jrr.best_matchedI{
        _,OK :=  alike_index.seqs[libseqnum]
        if 0==libseqnum||OK{continue}
        start := jrr.best_startsI[libseqnum]
        mm := jrr.best_mmI[libseqnum]
        add_one_aln(ofp, libseqnum,start,n_matched,mm)
        nwt ++
    }

    for libseqnum, n_matched := range jrr.best_matchedM{
        start := jrr.best_startsM[libseqnum]
        mm := jrr.best_mmM[libseqnum]
        add_one_aln(ofp, libseqnum,start,n_matched,mm)
        nwt ++
    }
    tp1 := time.Now()

    if verbose {println("RUN-Nread", run_no," TP_chan", tp0.Sub(tps), "  TP_write", tp1.Sub(tp0))}
    my_worker++
  }
  fmt.Println("Finished Alignment of Read", util.IfElseInt(is_R2,2,1) ,"FASTQ.gz File.")
}

func match_simple(q,r string,st,minM int)( ma,mm int ){
  match_Rstart := 0
  match_Rend := len(r)
  if st > 0{match_Rstart = st}
  if st+len(q) < len(r) { match_Rend = st+len(q) }
  maxMM := match_Rend - match_Rstart - minM

  ma=0
  if maxMM < 0 {
    mm=999999
    return
  }
  mm=0
  for Rii := match_Rstart; Rii < match_Rend; Rii++{
    qii :=0
    if st <0{ qii = Rii - st } else {qii = Rii - match_Rstart}

    if q[qii] == r[Rii] {ma++}else{
       mm++
       if mm > maxMM{
          ma=0
          mm=999999
          break
       }
    }
  }
  return
}


func worker(jobin chan jobdef, jobout chan jobres, refseq1 string,max_lib_no int, alike_index alike_seq_index, novalid_chrs int, no_chrs []bool, no_seqs_for_seqs []int){
    for{
        job, hasmore := <-jobin
        if ! hasmore{break}
        one_read_match(refseq1,max_lib_no, alike_index, job, novalid_chrs, no_chrs, no_seqs_for_seqs, jobout)
    }
}

// assumption: all of the ref sequences have the same length.
func alike_match_by_ints(q, refseq1 string, alike_index alike_seq_index, max_lib_no, novalid_chrs int, no_chrs []bool) (best_match_start_at_chros, best_matched_no, best_mismatched_no []int){
  len_ref := len(alike_index.ints[0][0])
  best_match_start_at_chros = make([]int, max_lib_no+1)
  best_matched_no = make([]int, max_lib_no+1)
  best_mismatched_no = make([]int, max_lib_no+1)
  for ii := range best_mismatched_no{best_mismatched_no[ii]=9999999}
  cur_best_scores := make([]int, max_lib_no+1)
  pilot_run_try_start_in_ref := -9999
  pilot_run_best_score_seq1 := 0
  min_qlen_for_run:= -9999
  for realrun := 0; realrun < 2; realrun ++{
    for try_start_in_refX :=-len(q); try_start_in_refX < len_ref; try_start_in_refX++{
      is_pilot_run := false
      try_start_in_ref := try_start_in_refX
      if try_start_in_refX == -len(q){
         if pilot_run_try_start_in_ref < -1000{continue}
         try_start_in_ref = pilot_run_try_start_in_ref
         is_pilot_run = true
      }else if realrun==1 && try_start_in_refX == pilot_run_try_start_in_ref{continue}

      rs := 0
      if try_start_in_ref >=0{ rs = try_start_in_ref }
      re := try_start_in_ref+len(q)
      if re > len_ref{re = len_ref}
      if re - rs < min_qlen_for_run{ continue }

      qs := 0
      if try_start_in_ref <0 {qs = -try_start_in_ref}
      qe := len(q)
      if qe - qs > re - rs { qe = qs + re - rs }

      q_sub := q[qs:qe]
      if realrun == 0{
         my_score :=0
         for qidx, qbase := range q_sub{
           refbase := rune(refseq1[qidx+rs])
           if refbase == qbase{ my_score ++}
         }
         if my_score > pilot_run_best_score_seq1{
           pilot_run_best_score_seq1 = my_score
           pilot_run_try_start_in_ref = try_start_in_ref
         }
         continue
      }

      chro_Match_counts := make([]int, max_lib_no+1,max_lib_no+1)
      chro_MisMatch_counts := make([]int, max_lib_no+1,max_lib_no+1)
      done_in_MisMatch_mode := 0
      stopped_try := make([]bool,  max_lib_no+1,max_lib_no+1)
      bad_gone := novalid_chrs

      var qbaseint int
      var stopped_in_middle bool
      for qidx, qbase := range q_sub{
        ridx := rs + qidx
        remain_len := len(q_sub)-qidx

        if qbase == 'A'{
          qbaseint = 0
        } else if qbase == 'C'{
          qbaseint = 1
        } else if qbase == 'G'{
          qbaseint = 2
        } else if qbase == 'T'{
          qbaseint = 3
        } else if qbase == 'N'{
          qbaseint = -888
        }

        if qbaseint>=0{
          matching_idx := alike_index.ints[0][qbaseint][ridx]
          if len(matching_idx) > max_lib_no/2{
            done_in_MisMatch_mode ++
            mismatching_idx := alike_index.ints[1][qbaseint][ridx]
            for _, toadd_chr := range mismatching_idx{
              if stopped_try[toadd_chr] {continue}
              newmisma := chro_MisMatch_counts[toadd_chr]
              vote_from_match := chro_Match_counts[toadd_chr]
              if (done_in_MisMatch_mode - newmisma) + vote_from_match + remain_len < cur_best_scores[toadd_chr]{
                stopped_try [toadd_chr] = true
                bad_gone ++
              }else{
                chro_MisMatch_counts[toadd_chr] = newmisma+1
              }
            }
          }else{
            for _, toadd_chr := range matching_idx{
              if stopped_try[toadd_chr] {continue}
              newvote := chro_Match_counts[toadd_chr]
              vote_from_misma := done_in_MisMatch_mode - chro_MisMatch_counts[toadd_chr]
              if newvote + vote_from_misma + remain_len < cur_best_scores[toadd_chr]{
                stopped_try [toadd_chr] = true
                bad_gone ++
              }else{
                chro_Match_counts[toadd_chr] = newvote+1
              }
            }
          }
        }
        if bad_gone >max_lib_no {
           stopped_in_middle =true
           break
        }
      }
      if stopped_in_middle{continue}

      // Must be real run. Can be pilot or not.
      for chro_i, chro_matched := range chro_Match_counts{
        if stopped_try[chro_i] || no_chrs[chro_i] { continue }
        chro_matched += (done_in_MisMatch_mode - chro_MisMatch_counts[chro_i])
        chro_misma := len(q_sub) - chro_matched
        if chro_matched > best_matched_no[chro_i] || (chro_matched == best_matched_no[chro_i] && chro_misma < best_mismatched_no[chro_i]) {
          best_matched_no[chro_i] = chro_matched
          best_match_start_at_chros[chro_i] = try_start_in_ref
          best_mismatched_no[chro_i] = chro_misma
        }
      }

      if is_pilot_run {
        cur_best_scores = chro_Match_counts
        for chr_no, misma := range chro_MisMatch_counts{ cur_best_scores[chr_no] += ( done_in_MisMatch_mode - misma ) }
        min_qlen := 88888
        for qchr, qv := range cur_best_scores{
          if qv < min_qlen && !no_chrs[qchr] {min_qlen = qv}
        }
        min_qlen_for_run = min_qlen
      }
    }
  }
  return
}




// ref sequences can have diff lengths.
// ref sequences are longer than query
func alike_match_by_seqs(q string, alike_index alike_seq_index, max_lib_no int, seq_nos[]int) (best_match_start_at_chros, best_matched_no, best_mismatched_no map[int]int){
  best_match_start_at_chros = make(map[int]int)
  best_matched_no = make(map[int]int)
  best_mismatched_no = make(map[int]int)

  for _,Rno := range seq_nos{
    Rseq := alike_index.seqs[Rno]
    best_mismatched := 999999
    best_matched := 0
    best_match_start_at := -999999

    for this_start:= - len(q)+1; this_start < len(Rseq) - 1; this_start++{
      my_match, my_misma := match_simple(q, Rseq, this_start, best_matched)
      if my_match > best_matched || (my_match == best_matched && my_misma < best_mismatched) {
         best_matched = my_match
         best_mismatched = my_misma
         best_match_start_at = this_start
      }
    }
    best_match_start_at_chros[Rno] = best_match_start_at
    best_matched_no[Rno] = best_matched
    best_mismatched_no[Rno] = best_mismatched
  }
  return
}

func add_one_rec(ofp  *bufio.Writer, q string){
  qbyte := []byte(q)
  if len(qbyte) <max_read_len{
    padding := bytes.Repeat([]byte{0x20}, max_read_len-len(qbyte))
    qbyte = append(qbyte, padding...)
  }
  ofp.Write(qbyte)
}

func add_one_aln(ofp  *bufio.Writer, libseqnum,start,n_matched,mm int){
  ii := int16(libseqnum)
  binary.Write(ofp, binary.LittleEndian, &ii)
  ii = int16(start)
  binary.Write(ofp, binary.LittleEndian, &ii)
  ii = int16(n_matched)
  binary.Write(ofp, binary.LittleEndian, &ii)
  ii = int16(mm)
  binary.Write(ofp, binary.LittleEndian, &ii)
}

type jobres struct{
    best_startsI,best_matchedI, best_mmI []int
    best_startsM, best_matchedM, best_mmM map[int]int
}

type jobdef string

func one_read_match(refseq1 string, max_lib_no int, alike_index alike_seq_index, qseq jobdef, novalid_chrs int, no_chrs[]bool, seq_nos []int, out chan jobres){
  best_startsI, best_matchedI, best_mmI := alike_match_by_ints(string(qseq), refseq1, alike_index, max_lib_no, novalid_chrs, no_chrs)
  best_startsM, best_matchedM, best_mmM := alike_match_by_seqs(string(qseq), alike_index, max_lib_no, seq_nos)
  out <- jobres{  best_startsI, best_matchedI, best_mmI , best_startsM, best_matchedM, best_mmM }
}
