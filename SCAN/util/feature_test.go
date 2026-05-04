package util

import "testing"
import "fmt"

func TestFeatures(t *testing.T){
	println("loading annotation...")
	include_introns := false
	feature_idx, err := LoadFeatureTable("/export/share/elvis/bioinf/liao/index/annotations/hg19-RefSeq.txt", include_introns)
	println("finished loading")
	if err != nil { panic(err) }

	hit_features := feature_idx.Lookup("chr1", 14362, 18366)
	for _,h :=range hit_features{
		fmt.Println(h.Key())
	}

	println ("Gene Neighbour between 100287654,100133331:", feature_idx.IsGeneNeighbours("100287654","100133331"))
	println ("Gene Neighbour between 148398,26155:", feature_idx.IsGeneNeighbours("26155","148398"))
	println ("Gene Neighbour between 100287654,26155:", feature_idx.IsGeneNeighbours("26155","100287654"))

	println("finished Features")
	println()
	println()

	features, _, err:= LoadFeatures("test-shuffle-features.SAF2",false)
	if err != nil { panic(err) }

	distinct_genes := make(map[string][]*GeneFeature)
	for _, ff := range features{
		_, OK := distinct_genes[ff.GeneID]
		if !OK{ distinct_genes[ff.GeneID]=make([]*GeneFeature,0) }
		distinct_genes[ff.GeneID]=append(distinct_genes[ff.GeneID], ff)
	}

	for _, fli := range distinct_genes{
		fli = MergeSortFeatures(fli)
		for fei, fe := range fli{
			if (fe.GeneID == "g1" && fei ==0 && fe.Stop != 1800) ||
			   (fe.GeneID == "g1" && fei ==1 && fe.Stop != 2200) ||
			   (fe.GeneID == "g2" && fei ==2 && fe.Start != 1801) {
				panic("MergeSort result is unright.")
			}
			println(fe.GeneID, " = ",fe.Chro , " : ", fe.Start, " ~ ", fe.Stop)
		}
	}
}
