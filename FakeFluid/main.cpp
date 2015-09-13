#include <functional>
#include <iostream>
#include <cstdint>

#include "../akari3/hlib/vec3.h"
#include "../akari3/hlib/image.h"

#include "perlin.h"

const int NX = 128;
const int NY = 128;
const int NZ = 128;

static float field[NX * NY * NZ];
static float field2[NX * NY * NZ];

using namespace hstd;

static akari3::PerlinImprovedNoise perlin;

double P1(const Float3& pos) {
    return perlin.noise(pos.x / 16.0f, pos.y / 16.0f, pos.z / 16.0f);
}
double P2(const Float3& pos) {
    return perlin.noise(pos.x / 16.0f + 100, pos.y / 16.0f + 200, pos.z / 16.0f + 300);
}
double P3(const Float3& pos) {
    return perlin.noise(pos.x / 16.0f - 100, pos.y / 16.0f - 200, pos.z / 16.0f - 300);
}

enum Dimension {
    DimX = 0,
    DimY = 1,
    DimZ = 2
};

template<unsigned Dim>
double dPdK(std::function<double(const Float3&)> f, const Float3& pos) {
    const double eps = 1e-5;
    Float3 pos0 = pos;
    Float3 pos1 = pos;

    const int idx = Dim;

    pos0[idx] -= eps;
    pos1[idx] += eps;

    return (f(pos1) - f(pos0)) / (2 * eps);
}

Float3 velocity(const Float3& pos) {
    const double dP3dY = dPdK<DimY>(P3, pos);
    const double dP2dZ = dPdK<DimZ>(P2, pos);
    const double dP1dZ = dPdK<DimZ>(P1, pos);
    const double dP3dX = dPdK<DimX>(P3, pos);
    const double dP2dX = dPdK<DimX>(P2, pos);
    const double dP1dY = dPdK<DimY>(P1, pos);

    return Float3(dP3dY - dP2dZ, dP1dZ - dP3dX, dP2dX - dP1dY);
}

unsigned int I(unsigned int ix, unsigned int iy, unsigned int iz) {
    return iz * (NX * NY) + iy * NX + ix;
}

void init() {
    for (int iz = 0; iz < NZ; ++iz) {
        for (int iy = 0; iy < NY; ++iy) {
            for (int ix = 0; ix < NX; ++ix) {
                const float x = (float)(ix + 0.5f) / 128.0;
                const float y = (float)(iy + 0.5f) / 128.0;
                const float z = (float)(iz + 0.5f) / 128.0;

                const float l = (Float3(x, y, z) - Float3(0.5f, 0.5f, 0.5f)).length();

                const bool b = (((ix / 8)%2 + (iy/8)%2 + (iz/8)%2) % 2 == 0) && l < 0.35f;

                field[I(ix, iy, iz)] = b;
            }
        }
    }
    /*
    for (int iz = 0; iz < NZ; ++iz) {
        for (int iy = 32; iy < NY - 32; ++iy) {
            for (int ix = 0; ix < NX; ++ix) {
                const float x = (float)ix / 256.0;
                const float y = (float)iy / 256.0;
                const float z = (float)iz / 256.0;

                const bool b = ((ix / 16)%2 + (iy/16)%2 + (iz/16)%2) % 2 == 0;

                field[I(ix, iy, iz)] = b;
            }
        }
    }
    */
}

void clear(float *f) {
    for (int iz = 0; iz < NZ; ++iz) {
        for (int iy = 0; iy < NY; ++iy) {
            for (int ix = 0; ix < NX; ++ix) {
                f[I(ix, iy, iz)] = 0;
            }
        }
    }
}

Float3 is2ws(const Int3& is_pos) {
    return (is_pos + Float3(0.5f, 0.5f, 0.5f));
}

Float3 uvw(const Float3& ws_pos) {
    // world space -> uvw space
    float u = ws_pos.x / (float)NX;
    float v = ws_pos.y / (float)NY;
    float w = ws_pos.z / (float)NZ;

    // texture wraping
    // repeat
    u = u - floor(u);
    v = v - floor(v);
    w = w - floor(w);

    return Float3(u, v, w);
}


float sample(float* field, const Float3& ws_pos) {
    const Float3 uv_pos = uvw(ws_pos);
    const Float3 is_pos = times(uv_pos, Float3(NX, NY, NZ));

    float ix = is_pos.x;
    float iy = is_pos.y;
    float iz = is_pos.z;
    
    int ix1 = (int)(ix + 0.5f) % NX;
    int ix0 = (ix1 - 1); if (ix0 < 0) ix0 += NX;
    float s0 = ix - (ix0 + 0.5f); if (s0 < 0) s0 += NX;
    float s1 = 1 - s0;

    int iy1 = (int)(iy + 0.5f) % NY;
    int iy0 = (iy1 - 1); if (iy0 < 0) iy0 += NY;
    float t0 = iy - (iy0 + 0.5f); if (t0 < 0) t0 += NY;
    float t1 = 1 - t0;

    int iz1 = (int)(iz + 0.5f) % NZ;
    int iz0 = (iz1 - 1); if (iz0 < 0) iz0 += NZ;
    float u0 = iz - (iz0 + 0.5f); if (u0 < 0) u0 += NZ;
    float u1 = 1 - u0;

    return
        s0 * t0 * u0 * field[I(ix0, iy0, iz0)] + 
        s0 * t0 * u1 * field[I(ix0, iy0, iz1)] + 
        s0 * t1 * u0 * field[I(ix0, iy1, iz0)] + 
        s0 * t1 * u1 * field[I(ix0, iy1, iz1)] + 
        s1 * t0 * u0 * field[I(ix1, iy0, iz0)] + 
        s1 * t0 * u1 * field[I(ix1, iy0, iz1)] + 
        s1 * t1 * u0 * field[I(ix1, iy1, iz0)] + 
        s1 * t1 * u1 * field[I(ix1, iy1, iz1)];
}

void advect(float* field, float* next_field, float dt) {
    std::cout << "advect" << std::endl;
    
#pragma omp parallel for schedule(dynamic, 1)
    for (int iz = 0; iz < NZ; ++iz) {
        for (int iy = 0; iy < NY; ++iy) {
            for (int ix = 0; ix < NX; ++ix) {
                const Float3 ws_pos = is2ws(Int3(ix, iy, iz));
                const Float3 vel = velocity(ws_pos * 32) / 2.0f + velocity(ws_pos * 2) /* + Float3(0.05, 0, 0) */;
                const Float3 prev_ws_pos = ws_pos - vel * dt;
                const float prev_value = sample(field, prev_ws_pos);

                // std::cout << vel << " ";

                next_field[I(ix, iy, iz)] = prev_value;
            }
        }
    }

}

struct Header {
    uint32_t magic;
    uint32_t X;
    uint32_t Y;
    uint32_t Z;

    Header(int x, int y, int z) : magic(0x01), X(x), Y(y), Z(z) {
    }
};

void save(float* f, const char* filename) {
    FILE* fp = fopen(filename, "wb");
    if (fp == NULL) {
        std::cout << "Error: save()" << std::endl;
        return;
    }

    Header header(NX, NY, NZ);
    fwrite(&header, sizeof(Header), 1, fp);

    fwrite(f, sizeof(float), NX * NY * NZ, fp);

    fclose(fp);
}

int main() {
    init();

    float* current = field2;
    float* prev = field;

    for (int i = 0; i < 5; ++i) {
        clear(current);
        advect(prev, current, 10.0f);
        float* tmp = current;
        current = prev;
        prev = tmp;
    }

    save(prev, "volume.float");
    return 0;
}