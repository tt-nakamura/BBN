#ifndef __nr_h__
#define __nr_h__

#include<cstdlib>
#include<string>
#include<iostream>

#include "Vec.h"
#include "Mat.h"

typedef double DP;

template<class T>
inline const T SIGN(const T &a, const T &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

template<class T>
inline const T MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}

template<class T>
inline const T MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}

template<class T>
inline const T SQR(const T a) {return a*a;}

inline void nrerror(const char *s)
{ std::cerr << s << std::endl; exit(1); }

#endif // __nr_h__
