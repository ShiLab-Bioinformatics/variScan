package util

type MyError struct{
	msg string
}

func NewError(msg string) *MyError{
	ret := new(MyError)
	ret.msg=msg
	return ret
}

func (self * MyError) Error() string{
	return self.msg
}
