package util

import "strings"
import "fmt"
import "math/rand"
import "errors"

func ContainStr(sli []string, target string) bool {
	for _, os := range sli {
		if os == target {
			return true
		}
	}

	return false
}

func MatchPrefix(onestr string, prefix_pool []string) (prefix string, found bool) {
	for _, px := range prefix_pool {
		if strings.HasPrefix(onestr, px) {
			found = true
			prefix = px
			return
		}
	}
	found = false
	return
}

const hash_prime = 2806197313
const hash_prime2 = 1644611641

func HashCode(str string) uint {
	ret := uint(0)
	for svi, sv := range str {
		ret ^= uint(sv) * hash_prime
		ret ^= uint(svi+1) * hash_prime2
		ret = ((ret & 0xffff) << 16) ^ ((ret & 0xffff0000) >> 16)
	}

	return ret
}

func JoinInt(split string, ints []int) string{
	ret:=""
	for r1, st := range ints{
		ret += fmt.Sprintf("%d",st)
		if r1 < len(ints)-1{ ret += split }
	}

	return ret
}

func JoinStr(split string, strs []string) string{
	ret:=""
	for r1, st := range strs{
		ret += st
		if r1 < len(strs)-1{ ret += split }
	}

	return ret
}

func MinStr(a, b string) string{
	if a<b{return a}
	return b
}

func MaxStr(a, b string) string{
	if a>b{return a}
	return b
}

var letterRunes = []rune("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890")
func RandomString(length int) string{
    b := make([]rune, length)
    for i := range b {
        b[i] = letterRunes[rand.Intn(len(letterRunes))]
    }
    return string(b)
}


var baseRunes = []rune("ATGCN")
func RandomBases(length int) string{
    b := make([]rune, length)
    for i := range b {
        b[i] = baseRunes[rand.Intn(len(baseRunes))]
    }
    return string(b)
}

var base4Runes = []rune("ATGC")
func Random4Bases(length int) string{
    b := make([]rune, length)
    for i := range b {
        b[i] = base4Runes[rand.Intn(len(base4Runes))]
    }
    return string(b)
}

func Map2Slice(m map[string] bool) []string{
	ret := make([]string,0)
	for s := range m{
		ret=append(ret,s)
	}
	return ret
}

func Phred64_to_33(s string) string{
	ret := ""
	for _, c := range s{
		if c - 64 < 1{
			panic(fmt.Sprintf("Char %d is not a phred-64 score!", c))
		}
		ret +=  string(rune(int(c) - (64 - 33)))
	}
	return ret
}

func HammingDistance(a,b string)(dist int, err error){
  if len(a) != len(b) { return -1, errors.New("Input strings must have same length.") }
  distance := 0
  for i:= range a{
    if a[i] != b[i] {
       distance++
    }
  }
  return distance, nil
}

