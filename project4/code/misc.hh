#ifndef MISC_HH
#define MISC_HH

#include <cfloat>

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

typedef long long int integer;
typedef long double number;

#endif // MISC_HH
