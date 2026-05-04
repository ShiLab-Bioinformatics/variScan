package util

import "fmt"

func AsciiColorStr(clr int) string{
	return fmt.Sprintf("%c[3%dm", 0x1b, clr)
}
