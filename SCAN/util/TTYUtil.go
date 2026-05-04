package util

/*
extern void set_pty_canon(int is_canon, int echo);
*/
import "C"
import "os"
import "fmt"
import "time"

func SetTtyCanonical(is_canon, has_echo bool) {
	is_canon_i := 0
	if is_canon {
		is_canon_i = 1
	}

	has_echo_i := 0
	if has_echo {
		has_echo_i = 1
	}

	C.set_pty_canon(C.int(is_canon_i), C.int(has_echo_i))
}


const CODE_STAT_MOUSE_EVENT = 1

func consoleCodeMachine(nch byte, stat * int, inbuff *[]byte, last_time *time.Time){
	fmt.Printf("RECV : %d (%c), Time_elapse=%d\n", nch, nch, time.Now().UnixNano() - last_time.UnixNano())
	if (*stat) == CODE_STAT_MOUSE_EVENT {
		(* inbuff) = append(*inbuff, nch)
		if len(*inbuff) == 5{
			(*stat) = 0
			fmt.Printf("COORDINATE=% 3d,% 3d, BUT=%X\n", (*inbuff)[3] - 32, (*inbuff)[4] - 32, (*inbuff)[2]);
		}
	}else if nch == 033 {
		// ESC
		(*stat) = CODE_STAT_MOUSE_EVENT
		(*inbuff) = make([]byte, 0)
		(*last_time) = time.Now()
	}else{
	//	fmt.Printf("RECV : %d\n", nch)
	}
}

func TestTty() {
	var stat int
	var consoleBuf []byte
	var lastTime time.Time
	for {
		inBuf := make([]byte, 1)
		rlen, err := os.Stdin.Read(inBuf)
		if rlen > 0 {
			consoleCodeMachine(inBuf[0], &stat, &consoleBuf, &lastTime)
		}
		if err != nil {
			println(err)
			break
		}
	}
}
