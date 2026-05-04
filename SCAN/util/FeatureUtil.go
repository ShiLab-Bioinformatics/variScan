package util

import "fmt"
import "sort"
import "strconv"
import "strings"

type GeneFeature struct{
	Chro string
	Start , Stop int // Stop is inclusive (len = Stop - Start + 1)
	GeneID string
	IsNegativeStrand bool
	ExonID int
	FeatureIndex int
}

type FeatureIndex struct{
	position_table map[string]map[int][]*GeneFeature
	neighbouring_gene_table map[string]bool	// tab[small_gene_name][large_gene_name] = 1
	bucket_size int
	RawFeatures []*GeneFeature
}

type FeatureReadMatching struct{
	Feature GeneFeature
	ChroStart, ChroStop, ReadStart, ReadLen int // ChroStop is inclusive (len = Stop - Start + 1)
}

type SortingFeature []*GeneFeature

func (s SortingFeature) Swap(i,j int){
	s[i], s[j] = s[j], s[i]
}

func (s SortingFeature) Less(i, j int) bool{
	if s[i].Start < s[j].Start {return true}
	return false
}

func (s SortingFeature) Len() int{
	return len(s)
}

func MergeSortFeatures(gfs []*GeneFeature) []*GeneFeature{
	var ret SortingFeature
	ret = gfs
	sort.Sort(ret)
	retm := make([]*GeneFeature,0)
	for ri, rnew := range ret{
		if ri >0{
			rlast := retm[len(retm)-1]
			if rnew.Start <= rlast.Stop && rnew.Stop > rlast.Stop{
				rlast.Stop = rnew.Stop
			} else if rnew.Start > rlast.Stop{
				retm =  append(retm,  rnew)
			}
		}else{
			retm = append(retm, rnew)
		}
	}

	return retm
}

type SortedFeatureReadMatching []FeatureReadMatching

func (mas SortedFeatureReadMatching) Len() int{
	return len(mas)
}

func (mas SortedFeatureReadMatching) Less(i,j int) bool{
	return mas[i].ChroStart < mas[j].ChroStart
}

func (mas SortedFeatureReadMatching) Swap(i,j int) {
	tm := mas[i]
	mas[i]=mas[j]
	mas[j]=tm
}

func UniqueFeatureLength(mas []FeatureReadMatching) int{
	if len(mas)<1 {return 0}
	ret := 0

	/*
	sort.Sort(SortedFeatureReadMatching(mas))

	gapstack = append(gapstack, mas[0])
	for mai, ma := range mas{
		if mai < 1{continue}
		top_last := gapstack[ len(gapstack) - 1 ]
		if top_last.ReadStart + top_last.ReadLen < ma.ReadStart {
			gapstack = append(gapstack, ma)
		} else if top_last.ReadStart + top_last.ReadLen < ma.ReadStart + ma.ReadLen {
			top_last.ReadLen = ma.ReadStart + ma.ReadLen - top_last.ReadStart
			gapstack[ len(gapstack) - 1 ] = top_last
		}
	}

	for _, ma := range gapstack{
		ret += ma.ReadLen
	}
	*/

	chro_mas := make(map[string] []FeatureReadMatching)
	for _, ma := range mas{
		chro := ma.Feature.Chro
		if len( chro_mas[chro] )<1{ chro_mas[chro] = make([]FeatureReadMatching,0) }
		chro_mas[chro] = append(chro_mas[chro] , ma)
		//println(ma.ChroStop , ma.ChroStart, chro,  len( chro_mas[chro] ))
	}

	for _, this_mas := range chro_mas{
		sort.Sort(SortedFeatureReadMatching(this_mas))
		gapstack := make([]FeatureReadMatching, 0)
		gapstack = append(gapstack, this_mas[0])
		for mai, ma := range this_mas{
			if mai < 1{continue}
			top_last := gapstack[ len(gapstack) - 1 ]
			if top_last.ChroStop + 1 < ma.ChroStart {
				gapstack = append(gapstack, ma)
			} else if top_last.ChroStop < ma.ChroStop {
				top_last.ChroStop = ma.ChroStop
				gapstack[ len(gapstack) - 1 ] = top_last
			}
		}
		for _, ma := range gapstack{
			//println(ma.ChroStop , ma.ChroStart)
			ret += ma.ChroStop - ma.ChroStart + 1
		}
	}

	return ret
}

type Exon struct{
	Chro string
	Start,Stop int //Both are 1-base coordinates. Both inclusive in this exon. The same as in GTF format.
	IsNegative bool
}

type Transcript struct{
	GeneID, TranscriptID string
	Exons []Exon
}

// lookup_table: chro:pos1 <TAB> pos2.   pos1 and pos2 are both 1-base coordinates, both included in the exons. pos1 < pos2.
// value: 1 if positive; 2 if negative; bitwise-or 
func LoadJuncsInTranscripts(GTF_fname string, gtf_gene_id, gtf_transcript_id string)(lookup_table map[string]int, err error){
	txns, terr := LoadTranscripts(GTF_fname, gtf_gene_id, gtf_transcript_id)
	if terr!=nil{
		err = terr
		return
	}

	lookup_table = make(map[string]int)

	for _, tx := range txns{
		var exx Exon
		for exi, ex := range tx.Exons{
			if exi==0{
				exx = ex
				continue
			}

			if ex.Chro != exx.Chro || ex.IsNegative != exx.IsNegative {panic("Transcript jumpped!")}
			var left_val,right_val int
			if ex.IsNegative{
				right_val = exx.Start
				left_val = ex.Stop
			}else{
				left_val = exx.Stop
				right_val = ex.Start
			}
			if left_val >= right_val{panic(fmt.Sprintf("Order of exons is wrong: %d >= %d",  left_val, right_val))}
			look_key := ex.Chro+":"+strconv.Itoa(left_val)+"\t"+strconv.Itoa(right_val)
			OKv := lookup_table[look_key]
			if ex.IsNegative {lookup_table[look_key] = OKv|2 } else{lookup_table[look_key] = OKv|1 }
			exx = ex
		}
	}
	err = nil
	return
}

func LoadTranscripts(GTF_fname string, gtf_gene_id, gtf_transcript_id string) (txns []Transcript, err error){
	gfp,ferr := Qopen(GTF_fname)
	if ferr != nil{
		err=ferr
		return
	}
	defer gfp.Close()

	tns_table := make(map[string][]Exon) // key: GeneID+"@^@"+TranscriptID
	for{
		fli, aerr := gfp.Array()
		if aerr != nil{break}

		if len(fli)<9{continue}
		if fli[0][0]=='#'{continue}
		if fli[2]!="exon"{continue}

		chro := fli[0]
		exon_start,_ := strconv.Atoi(fli[3]) // Genomic start of the feature (inclusive), with a 1-base offset.
		exon_end,_ := strconv.Atoi(fli[4]) // Genomic end of the feature (inclusive), with a 1-base offset.
		is_neg := fli[6]=="-"
		if (fli[6]=="+") == is_neg {panic("The strand isn't + nor - !")}

		gene_name := strings.Split(fli[8], gtf_gene_id+" \"")[1]
		gene_name = strings.Split(gene_name, "\"")[0]

		transcript_name := strings.Split(fli[8], gtf_transcript_id+" \"")[1]
		transcript_name = strings.Split(transcript_name, "\"")[0]

		gene_txn_key := gene_name +"@^@"+transcript_name
		exon_list, OK := tns_table[gene_txn_key]
		if !OK{ exon_list = make([]Exon,0) }
		exon_list = append(exon_list, Exon{ chro, exon_start, exon_end, is_neg })
		tns_table[gene_txn_key] = exon_list
	}

	txns = make([]Transcript,0)
	for ky, exons := range tns_table{
		keyinf := strings.Split(ky,"@^@")
		gene_name := keyinf[0]
		transcript_name := keyinf[1]
		txns = append(txns, Transcript{ gene_name, transcript_name, exons })
	}
	err = nil
	return
}

func LoadFeatures(fname string, include_intron bool) (features []*GeneFeature, neighbour_gene_table map[string]bool, err error){
	features, neighbour_gene_table, err = LoadFeaturesEx(fname, include_intron, "gene_id")
	return
}

func LoadFeaturesEx(fname string, include_intron bool, gtf_gene_id string) (features []*GeneFeature, neighbour_gene_table map[string]bool, err error){
	af, xerr := Qopen(fname)

	if xerr != nil{
		err = xerr
	}else{
		features = make([]*GeneFeature,0)
		neighbour_gene_table = make(map[string]bool)
		old_gene_id := ""
		old_chro := ""
		exon_n := 1

		for{
			qarr, yerr := af.Array()
			if yerr != nil{ break }

			if qarr[0][0]=='#'{ continue }
			if len(qarr)<4{ continue }
			if qarr[0] == "GeneID" || qarr[0] == "Geneid"{ continue }

			newfeature := new(GeneFeature)
			if(len(qarr)>6){
				FeatureMode := qarr[2]
				if FeatureMode != "exon"{continue}
				pos1, _ := strconv.Atoi(qarr[3])
				pos2, _ := strconv.Atoi(qarr[4])

				newfeature.Chro = qarr[0]
				newfeature.Start = pos1
				newfeature.Stop = pos2
				newfeature.IsNegativeStrand = qarr[6] == "-"
				newfeature.GeneID = strings.Split(strings.Split(qarr[8], gtf_gene_id +" \"")[1], "\"")[0]
			}else{

				pos1, _ := strconv.Atoi(qarr[2])
				pos2, _ := strconv.Atoi(qarr[3])

				newfeature.Chro = qarr[1]
				newfeature.Start = pos1
				newfeature.Stop = pos2
				if strings.Contains( qarr[0], "=" ){
					newfeature.GeneID = strings.Split(qarr[0], "=")[1]
				}else{
					newfeature.GeneID = qarr[0]
				}
				newfeature.IsNegativeStrand = qarr[4] == "-"
			}
			if  old_gene_id!= newfeature.GeneID{ exon_n = 1 }
			newfeature.ExonID = exon_n

			exon_n ++

			newfeature.FeatureIndex = len(features)
			features = append(features, newfeature)

			if old_chro != newfeature.Chro{
				old_chro = newfeature.Chro
				old_gene_id = ""
			}

			if old_gene_id!= newfeature.GeneID{
				if old_gene_id != ""{
					nbkey :=fmt.Sprintf("%s::%s", MinStr(old_gene_id, newfeature.GeneID), MaxStr(old_gene_id, newfeature.GeneID) )
					neighbour_gene_table[nbkey] = true
				}

				old_gene_id = newfeature.GeneID
			}
		}

		if include_intron{
			genes := make(map[string] *GeneFeature)
			for _, f := range features{
				_, ok:=genes[f.GeneID]
				if ok{
					if f.Chro != genes[f.GeneID].Chro{
						f.Chro = "MULTI-CHROMOSOME"
					}
					f.Start = MinInt(f.Start, genes[f.GeneID].Start)
					f.Stop = MaxInt(f.Stop, genes[f.GeneID].Stop)
				}
				genes[f.GeneID] = f
			}
			features = make([]*GeneFeature,0)
			for _, f := range genes{
				features = append(features, f)
			}
		}
	}

	af.Close()

	if len(features)<1 { err = NewError("No feature was found in "+fname) }

	return
}

func (self *GeneFeature) Key() string{
	return fmt.Sprintf("%s:%s:%d~%d", self.GeneID, self.Chro, self.Start, self.Stop)
}

func LoadFeatureTable(fname string, include_intron bool) ( FeatureIndex *FeatureIndex , err error){
	features, neighbour_table, xerr := LoadFeatures(fname, include_intron)


	if xerr != nil{ err = xerr }else{
		FeatureIndex = CreateFeatureIndex()
		FeatureIndex.neighbouring_gene_table = neighbour_table
		FeatureIndex.RawFeatures = features

		for _, fe := range features{
			FeatureIndex.Append(fe.Chro, fe.Start, fe.Stop, fe)
		}
	}

	return
}

func CreateFeatureIndex() *FeatureIndex{
	ret := new(FeatureIndex)
	ret.position_table = make(map[string]map[int][]*GeneFeature)
	ret.bucket_size = 10240
	return ret
}

func (self *FeatureIndex) Append(chro string, start, end int, feat * GeneFeature){
	_, ok := self.position_table[chro]
	if !ok{ self.position_table[chro] = make(map[int][]*GeneFeature) }

	bucket_0 := start - start % self.bucket_size
	bucket_n := end - end % self.bucket_size

	for xk1 := bucket_0; xk1 <= bucket_n; xk1+=self.bucket_size {
		_, ok = self.position_table[chro][xk1]
		if !ok{ self.position_table[chro][xk1] = make([]*GeneFeature,0) }

		self.position_table[chro][xk1] = append(self.position_table[chro][xk1], feat)
	}
}

func (self *FeatureIndex) Features(sec MappedSection) []*GeneFeature{
	return self.Lookup(sec.Chro, sec.Pos, sec.Pos + sec.ChromosomeLength - 1)
}

func (self *FeatureIndex) Lookup(chro string, start, stop int) []*GeneFeature{
	retmap := make(map[string]*GeneFeature)

	_, ok := self.position_table[chro]
	//fmt.Printf("FOUND = %t\n", ok);

	if ok{
		//fmt.Printf("FOUND_LEN = %d\n", len(tabs));

		bucket_0 := start - start % self.bucket_size
		bucket_n := stop - stop % self.bucket_size


		for xk1 := bucket_0; xk1 <= bucket_n; xk1+=self.bucket_size {
			arrv, ok := self.position_table[chro][xk1]
			if ok{
				for _, fep := range arrv{
					fky := fep.Key()
					retmap[fky] = fep
				}
			}
		}
	}

	ret := make([]*GeneFeature, 0)

	for _, fptr := range retmap{
		if start <= fptr.Stop && stop >= fptr.Start{
			ret = append(ret, fptr)
		}
	}
	return ret
}

func (self *FeatureIndex) IsNeighbourFeatures(f1, f2 GeneFeature) bool{
	return self.IsGeneNeighbours(f1.GeneID, f2.GeneID)
}
func (self *FeatureIndex) IsGeneNeighbours(id1, id2 string) bool{
	id_min := MinStr(id1, id2)
	id_max := MaxStr(id1, id2)

	return self.neighbouring_gene_table[fmt.Sprintf("%s::%s", id_min, id_max)]

}

func (self *FeatureIndex) GetGeneCount() int{
	return len(self.neighbouring_gene_table)
}
