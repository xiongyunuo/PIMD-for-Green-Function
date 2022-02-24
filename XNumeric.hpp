#ifndef XNUMERIC_HPP
#define XNUMERIC_HPP

#include <vector>
#include <cstdlib>
#include <iostream>

#define XNUMERRORCODE 10

#define XERROR(str) {\
std::cerr << str << std::endl;\
std::exit(XNUMERRORCODE);\
}

typedef double XNum;
typedef std::vector<XNum> XVecD;
typedef std::vector<std::vector<XNum> > XMatD;
typedef XNum XOneArguFunc(XNum);

#endif