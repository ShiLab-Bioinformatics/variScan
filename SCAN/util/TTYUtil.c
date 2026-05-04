#include <stdio.h>
#include <stdlib.h>
#include <termios.h>
#include <unistd.h>

void set_pty_canon(int is_canon, int echo){
	struct termios attr;

	tcgetattr(fileno(stdin), &attr);
	if(is_canon)
		attr.c_lflag |= ICANON;
	else
		attr.c_lflag &= ~ICANON;

	if(echo)
		attr.c_lflag |= ECHO;
	else
		attr.c_lflag &= ~ECHO;

	tcsetattr(fileno(stdin), TCSANOW, &attr);
}
