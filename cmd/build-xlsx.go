package main

import "github.com/xuri/excelize/v2"
import "github.com/ShiLab-Bioinformatics/VariantAlign/SCAN/util"
import "encoding/csv"
import "os"
import "fmt"
import "strings"
import "strconv"

func main() {
    qnames := make([]string,0)

    fileName := os.Args[1]
    outfile := os.Args[2]
    file, _ := os.Open(fileName)
    reader := csv.NewReader(file)

    for {
        record, err := reader.Read()
        if err != nil { break }

        if len(record) >= 2 {
            qnames = append(qnames,strings.ReplaceAll(record[1]," ","_"))
        }else{
            panic("Bad reference format")
        }
    }

    qf,_:= util.Qopen("STDIN://")
    counts :=make(map[string]int)
    for{
        fl,err := qf.Line()
        if err!=nil{break}
        fl=strings.TrimSpace(fl)
        if !strings.Contains(fl,"total_num_mapped_sequences 1 mapped_to"){continue} // we only need uniqly mapped read-pairs (i.e., mapped to 1 seq).
        qname := strings.Split(fl," ")[12] // yes it is there.
        counts[qname]++
    }

    f := excelize.NewFile()
    index, _ := f.NewSheet("Read-pair Counts")
    f.SetCellValue("Read-pair Counts", "A1", "Sequence name")
    f.SetCellValue("Read-pair Counts", "B1", "Number of read-pairs")
    for qi,qname := range qnames{
        f.SetCellValue("Read-pair Counts", fmt.Sprintf("A%d", qi+2), qname)
        f.SetCellValue("Read-pair Counts", fmt.Sprintf("B%d", qi+2), strconv.Itoa(counts[qname]))
    }
    f.SetActiveSheet(index)
    f.DeleteSheet("Sheet1")
    if err := f.SaveAs(outfile); err != nil {
        fmt.Println(err)
    }
}
