/*********************************************************************************
*
* This code transpiles OpenCL C to C++ in order to run it on CPU.
* 
* OpenCL vector initialization is completely different from C/C++.
* I was not able to make it compile with any kind of macros or C++ operator
* overloading. Even overloading the comma operator was tried.
*
*********************************************************************************/

#include <string>
#include <fstream>
#include <sstream>
#include <assert.h>

#include <iostream>

static bool verbose = false;
static bool very_verbose = false;
static const std::string VEC = "glm::vec";

/**
* Use preprocessor to include a vector lib compatible with OpenCL vector functions.
* Currently including GLM.
* Note swizzling can cause name collisions: https://glm.g-truc.net/0.9.1/api/a00002.html
*/
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

/**
* Define special OpenCL functions here.
* Note that some have no meaning in single thread (define dummies for these).
*/
void prependDefinitions(std::string& src) {
	std::string defs = "                                     \n\
		#include <stdio.h> // printf                         \n\
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
		float fract(float x, float* outFloor) {              \n\
			*outFloor = floor(x);                            \n\
			return x - (*outFloor);                          \n\
		}                                                    \n\
		using namespace glm;                                 \n\
		#define float3 vec3                                  \n\
		#define float2 vec2                                  \n\
	";
	src = defs + "\n\n" + src;
}

/**
*
*/
void skipNestedParanthesis(std::string& src, int* j, char* c, int* lineCount) {
	int k = (*j)+1;
	char d = src[k];
	while (d != ')') {
		if (d == '\n') (*lineCount)++;
		if (d == '(') skipNestedParanthesis(src, &k, &d, lineCount);
		k++;
		d = src[k];
	}
	if (very_verbose) std::cout << "Line " << *lineCount << ": Skip nested paranthesis " << src.substr(*j, k - (*j) + 1) << std::endl;
	*j = k;
	*c = src[*j];
}

/**
*
*/
void skipLineComment(std::string& src, int* j, char* c, int* lineCount) {
	int k = (*j)+1;
	char d = src[k];
	while (d != '\n' && d != '\0'/*EOF*/) {
		k++;
		d = src[k];
	}
	if (very_verbose) std::cout << "Line " << *lineCount << ": Skip line comment " << src.substr(*j-1, k - (*j) + 1) << std::endl;
	(*lineCount)++;
	*j = k;
	*c = src[*j];
}

/**
*
*/
void skipBlockComment(std::string& src, int* j, char* c, int* lineCount) {
	int k = (*j);
	char d = src[k];
	char e = src[k-1];
	int startLine = *lineCount;
	while (!(e == '*' && d == '/')) {
		if (d == '\n') (*lineCount)++;
		e = src[k];
		k++;
		d = src[k];
	}
	if (very_verbose) std::cout << "Line " << startLine << " to " << *lineCount << ": Skip block comment\n"
		<< src.substr(*j-1, k - (*j) + 2) << std::endl;
	*j = k;
	*c = src[*j];
}

/**
*
*/
bool validCSymbolCharacter(char a) {
	return ((a >= 48 && a <= 57/*numbers*/) 
		|| (a >= 65 && a <= 90/*uppercase*/) 
		|| (a >= 97 && a <= 122/*lowercase*/) 
		|| a == 95/*underscore*/);
}

/**
* Turn swizzles into function calls, which can be handled much better by GLM.
* Currently only .xy supported. PLEASE EXTENT!
*/
int replaceSwizzling(std::string& src, int i) {
	if (i > src.length()-3) return i;
	if (src[i]=='.' && src[i+1]=='x' && src[i+2]=='y') {
		src.insert(i+3, "()");
		return i+5; // continue after ()
	}
	return i;
}

/**
*
*/
void replaceVectorTypes(std::string& src) {

	// the cast like operator that comes before CL vector initializations helps when parsing
	const std::string clVec = "(floatX)";
	const std::string clVecCustom1 = "(RealX)";

	int lineCount = 1;
	char a = 0;
	char b = 0;
	int i = 0;
	b = src[i];
	while (b != '\0') {

		// Replace .xy by .xy() and set i and b to after ()
		b = src[i = replaceSwizzling(src, i)];

		// Get preceeding char, ignoring whitespace
		int lastNonWhitespaceIndex = -1;
		if (i > 0) {
			int j = i-1;
			a = src[j];
			while (a == ' ' || a == 9/*tab*/) {
				j--;
				if (j < 0) a = 0;
				else a = src[j];
			}
			lastNonWhitespaceIndex = j;
		}

		// Skip comments and count lines
		if (a == '/' && b == '/') skipLineComment(src, &i, &b, &lineCount);
		else if (a == '/' && b == '*') skipBlockComment(src, &i, &b, &lineCount);
		else if (b == '\n') lineCount++;

		// Detect opening paranthesis that is not part of function notation
		bool preceededByFuncName = validCSymbolCharacter(a);
		if (b == '(' && !preceededByFuncName) {
			int commaCount = 0;
			int j = i+1;
			char c = src[j];
			char d = src[j-1];

			// Search paranthesis contents for commas
			while (c != ')') {
				if (c == ',') commaCount++;
				else if (d == '/' && c == '/') skipLineComment(src, &j, &c, &lineCount);
				else if (d == '/' && c == '*') skipBlockComment(src, &j, &c, &lineCount);
				else if (c == '(') skipNestedParanthesis(src, &j, &c, &lineCount);

				// if skipNestedParanthesis executed at this point, j points to a closing ")"
				d = src[j];
				j++;
				c = src[j];

				// catch unclosed paranthesis
				if (j >= src.length()) {
					std::cerr << "CL source lacks corresponding closing paranthesis for line " << lineCount << std::endl;
					exit(1);
				}
			}

			if (commaCount >= 1) { // Found vector initialization
				//TODO this will probably also detect if(int i = 0, j = 0;;) as vec2. Prevent this!

				// Detect the cast like operator that appears before CL vector initializations
				// and extract vector size number as char
				char vecN = 0;
				if (lastNonWhitespaceIndex >= clVec.length()) {
					if ((src.substr(lastNonWhitespaceIndex - clVec.length() + 1, clVec.length()-2) == clVec.substr(0, clVec.length() - 2)
						&& src[lastNonWhitespaceIndex] == clVec[clVec.length()-1])
						||
						(src.substr(lastNonWhitespaceIndex - clVecCustom1.length() + 1, clVecCustom1.length()-2) == clVecCustom1.substr(0, clVecCustom1.length() - 2)
						&& src[lastNonWhitespaceIndex] == clVecCustom1[clVecCustom1.length()-1])) {
						vecN = src[lastNonWhitespaceIndex - 1];
					}
				}

				// (if cast like operator not found use comma count to determine vector size)
				if (verbose) std::cout << "Line " << lineCount << ": Found vec" << vecN ? vecN : (commaCount+1);
				if (verbose) std::cout << " ..." << a << src.substr(i,j-i) << c << "..." << std::endl;
				if (commaCount == 1 || commaCount == 2 || commaCount == 3) {

					// Insert constructor call before opening paranthesis (position i)
					const std::string r = VEC + (vecN ? vecN : (char)(48+commaCount+1));
					src.insert(i, r);

					// Continue at closing paranthesis (position j + length of inserted text)
					i = j + r.length() - 1;
				} else {

					std::cerr << "Transpiler cannot replace vector of size " << (commaCount+1) << "." << std::endl;
					std::cerr << "Detected here: " << " ..." << a << src.substr(i,j-i) << c << "..." << std::endl;
					exit(1);
				}
			}
		}
		i++;
		b = src[i];
	}
}

/**
*
*/
void replaceStructInitializers(std::string& src) {
	//TODO remove the cast-like operator
	// consider: typedef struct{ int a,b,c; } Dummy;
	// valid in C/C++:  return { 1,2,3 }
	// valid in OpenCL: return (Dummy){ 1,2,3 }
	// valid in both:   Dummy d = { 1,2,3 }
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
	std::ofstream out(std::string(args[1]) + ".transpiled.cpp");
	out << str;
	return 0;
}