
#include <R.h>
#include "die.h"

#define BUFSIZE 1024

void die(const char *format, ...)
{
    char buf[BUFSIZE];
    va_list arg;

    va_start(arg, format);
    vsnprintf(buf, BUFSIZE, format, arg);
    va_end(arg);
    buf[BUFSIZE - 1] = '\0';
    error(buf);
}

