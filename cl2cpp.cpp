#include <string>
#include <fstream>
#include <sstream>
#include <assert.h>

#include <iostream>

static bool verbose = false;
static bool very_verbose = false;
static const std::string VEC = "glm::vec";

void includeVectorLib(std::string& src) {
	std::string lib = "                                     \n\
		#define GLM_FORCE_SWIZZLE                           \n\
		#include \"glm/vec2.hpp\"                           \n\
		#include \"glm/vec3.hpp\"                           \n\
		#include \"glm/geometric.hpp\" // dot, cross...     \n\
		#include \"glm/trigonometric.hpp\" // sin, cos...   \n\
		#include \"glm/exponential.hpp\" // sqrt, log...    \n\
	";
	src = lib + "\n" + src;
}

void prependDefinitions(std::string& src) {
	std::string defs = "                                     \n\
		#include <fstream> // ofstream                       \n\
		#include <iomanip> // setw                           \n\
		#include <stdint.h> // uint32_t, uint64_t            \n\
		#define __constant                                   \n\
		#define __kernel                                     \n\
		#define __global                                     \n\
		#define uint uint32_t                                \n\
		#define ulong uint64_t                               \n\
		size_t get_global_id(uint dimindx) {                 \n\
			return 0;                                        \n\
		}                                                    \n\
		size_t get_global_size (uint dimindx) {              \n\
			return 1;                                        \n\
		}                                                    \n\
		unsigned int atomic_add(                             \n\
		volatile __global unsigned int *p ,                  \n\
		unsigned int val) {                                  \n\
			unsigned int old = *p;                           \n\
			*p += val;                                       \n\
			return old;                                      \n\
		}                                                    \n\
		#define float3 glm::vec3                             \n\
		#define dot glm::dot                                 \n\
		#define sign glm::sign                               \n\
		#define length glm::length                           \n\
		#define min glm::min                                 \n\
		#define xy xy()                                      \n\
	";
	src = defs + "\n\n" + src;
}

void skipNestedParanthesis(std::string& src, int* j, char* c, int lineCount) {
	int k = (*j)+1;
	char d = src[k];
	while (d != ')') {
		if (d == '(') skipNestedParanthesis(src, &k, &d, lineCount);
		k++;
		d = src[k];
	}
	if (very_verbose) std::cout << "skip " << src.substr(*j, k - (*j) + 1) << " at line " << lineCount;
	*j = k;
	*c = src[*j];
	if (very_verbose) std::cout << ", continue at " << src[*j] << src[*j+1] << src[*j+2] << "..." << std::endl;
}

void replaceVectorTypes(std::string& src) {
	int lineCount = 1;
	char a = 0;
	char b = 0;
	int i = 0;
	b = src[i];
	while (b != '\0') {
		if (i > 0) {
			int j = i-1;
			a = src[j];
			while (a == ' ' || a == 9/*tab*/) {
				j--;
				if (j < 0) a = 0;
				else a = src[j];
			}
		}
		bool funcName = ((a >= 48 && a <= 57/*numbers*/) 
			|| (a >= 65 && a <= 90/*uppercase*/) 
			|| (a >= 97 && a <= 122/*lowercase*/) 
			|| a == 95/*underscore*/);
		if (b == '(' && !funcName) {
			int commaCount = 0;
			int j = i+1;
			char c = src[j];
			while (c != ')') {
				if (c == '(') skipNestedParanthesis(src, &j, &c, lineCount);
				if (c == ',') commaCount++;
				j++;
				c = src[j];
			}
			if (commaCount >= 1) {
				if (verbose) std::cout << "found vec" << (commaCount+1) << " at line " << lineCount << ": ";
				if (verbose) std::cout << "..." << a << src.substr(i,j-i) << c << "..." << std::endl;
				if (commaCount == 1 || commaCount == 2 || commaCount == 3) {
					// This is a vector initialization
					const std::string r = VEC + std::to_string(commaCount+1);
					src.insert(i, r);
					i = j + r.length();
					continue;
				}
			}
		}
		if (b == '\n') lineCount++;
		i++;
		b = src[i];
	}
}

int main(int nargs, char* args[]) {
	assert(nargs>=2);
	std::ifstream in(args[1]);
	if (nargs>=3)
		if (std::string(args[2]) == "-v") verbose = true;
		else if (std::string(args[2]) == "-vv") very_verbose = verbose = true;
	std::stringbuf buf;
	in >> &buf;
	std::string str = buf.str();
	if (verbose) std::cout << "Read " << str.length() << " bytes" << std::endl;
	replaceVectorTypes(str);
	prependDefinitions(str);
	includeVectorLib(str);
	std::ofstream out(std::string(args[1]) + ".cpp");
	out << str;
	return 0;
}
