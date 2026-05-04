package util

import "math/rand"
import crand "crypto/rand"
import "fmt"
import "io"

func ChoiceStrExcept(strs []string, exc string) string{
	ret_idx := rand.Int() % (len(strs) - 1)
	for si :=0 ; si <= ret_idx; si++{
		if strs[si] == exc{
			return strs[ret_idx+1]
		}
	}
	return  strs[ret_idx]
}

func ChoiceStr(strs []string) string{
	rand_max := len(strs)
	rand_ind := rand.Int() % rand_max
	return strs[rand_ind]
}

func ChoiceInt(ints []int) int{
	rand_max := len(ints)
	rand_ind := rand.Int() % rand_max
	return ints[rand_ind]
}

func UUID() string {
	uuid := make([]byte, 16)
	_, _ = io.ReadFull(crand.Reader, uuid)
	// variant bits; see section 4.1.1
	uuid[8] = uuid[8]&^0xc0 | 0x80
	// version 4 (pseudo-random); see section 4.1.3
	uuid[6] = uuid[6]&^0xf0 | 0x40
	return fmt.Sprintf("%x-%x-%x-%x-%x", uuid[0:4], uuid[4:6], uuid[6:8], uuid[8:10], uuid[10:])
}


const b55str = "23456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnpqrstuvwxyz"
const b33str = "ABCDEFGHIJKLMNPQRSTUVWXYZ123456789"
func UXID() string{
	return rstr(b33str, 16)
}

func ULID() string{
	return rstr(b55str, 18)
}

func rstr(rb string, rl int) string{
	ret := ""
	for x:=0; x<rl; x++{
		v2 := make([]byte,3)
		crand.Reader.Read(v2)
		v22 := int(v2[1])*359+int(v2[0])+int(v2[2])*38281217
		ret += string ( rb[v22%len(rb)] )
	}
	return ret
}


