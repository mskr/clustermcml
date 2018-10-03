#ifndef LOG_H
#define LOG_H

#include <stdint.h> // uint32_t, uint64_t
#include <stdio.h> // printf, stdout

/**
* A logger like iostream
* but exception-free
*/
class Log {
public:
	static const char flush = '\0';
	Log& operator<<(char c) {
		if (c == flush) fflush(stdout);
		else printf("%c", c);
		return *this;
	}
	Log& operator<<(const char* str) {
		printf("%s", str);
		return *this;
	}
	Log& operator<<(int i) {
		printf("%d", i);
		return *this;
	}
	Log& operator<<(uint32_t u) {
		printf("%u", u);
		return *this;
	}
	Log& operator<<(uint64_t u) {
		printf("%llu", u);
		return *this;
	}
	Log& operator<<(double f) {
		printf("%f", f);
		return *this;
	}
};

static Log out;

#endif