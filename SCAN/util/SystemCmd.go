package util

import "os/exec"
import "io"

func System(cmd string) (return_code int, stdout string){
	out, err := exec.Command("/bin/bash", "-c", cmd).Output()
	return_code = 0
	if err != nil{ return_code = -1 }
	stdout = string(out)
	return
}

func ExecStream(cmdstr string) (return_code int, stdout io.ReadCloser){
	cmd := exec.Command("/bin/bash", "-c", cmdstr)
	return_code = 0

	var err error
	stdout, err = cmd.StdoutPipe()
	if err != nil{ return_code = -1 }
	err = cmd.Start()
	if err != nil{ return_code = -1 }
	return
}
