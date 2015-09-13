#ifndef _PERLIN_H_
#define _PERLIN_H_

#include <cmath>

#include "../akari3/hlib/random.h"

namespace akari3 {

/* imported from http://mrl.nyu.edu/~perlin/noise/ */

class PerlinImprovedNoise {
private:
public:
    PerlinImprovedNoise();
    double noise(double x, double y, double z);
};

};

#endif // _PERLINE_H_