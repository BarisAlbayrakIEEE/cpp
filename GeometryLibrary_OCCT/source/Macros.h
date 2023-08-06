#ifndef _Macros_HeaderFile
#define _Macros_HeaderFile

#include <cmath>
#include <string>
#include <array>
#include <vector>
#include <map>
#include <limits>

using arrayS2 = std::array<double, 2>;
using arrayS3 = std::array<double, 3>;
using arrayS4 = std::array<double, 4>;
using arrayS32 = std::array<std::array<double, 2>, 3>;
using arrayS33 = std::array<std::array<double, 3>, 3>;
using vectorInput1D = const std::vector<double, std::allocator<double>>&;
using vectorInput2D = const std::vector<std::vector<double, std::allocator<double>>>&;

#define ARGCOPY(className) const Handle(className)&
#define ARGMOVE(className) Handle(className)&&
#define OUTREF(className) Handle(className)&
#define OUTVAL(className) Handle(className)

#endif