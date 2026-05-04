package main

import "github.com/xuri/excelize/v2"
import "fmt"

func main() {
    f := excelize.NewFile()
    // Create a new sheet
    index, _ := f.NewSheet("Sheet1")
    // Set value of a cell
    f.SetCellValue("Sheet1", "A1", "Hello world.")
    // Set active sheet of the workbook
    f.SetActiveSheet(index)
    // Save spreadsheet by the given path
    if err := f.SaveAs("Book1.xlsx"); err != nil {
        fmt.Println(err)
    }
}
