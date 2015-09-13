#include <iostream>
#include <cmath>

#include "hlib\image.h"
#include "hlib\sampling.h"
#include "hlib\random.h"
#include "hlib\rt\ray.h"
#include "hlib\rt\ibl.h"

#include "hlib\memfile.h"
#include "hlib\hmath.h"

#include "hlib\rt\sphere.h"

#include "hlib\timer.h"
#include "hlib\rt\camera.h"
#include "hlib\rt\bbox.h"
#include "hlib\rt\objMesh.h"

#include "hlib\ext\bmpexporter.h"

using namespace hstd;

class SimpleVolume {
private:
    std::vector<float> density_;

    int Nx_, Ny_, Nz_;

    float max_density_;
    float max_scattering_;
    float max_transmittance_;
    float albedo_;

    Float3 org_; // cell(0, 0, 0)の端っこの点の世界座標系における座標
    Float3 scale_; // m / cell

    Float3 factor_;
    
    inline unsigned int I(unsigned int ix, unsigned int iy, unsigned int iz) const {
        return iz * (Nx_ * Ny_) + iy * Nx_ + ix;
    }
    
    inline Float3 is2ws(const Int3& is_pos) const {
        return (is_pos + Float3(0.5f, 0.5f, 0.5f));
    }

    inline Float3 uvw(const Float3& ws_pos) const {
        Float3 pos;
        
        pos = ws_pos - org_;

        Float3 u = times(factor_, pos);

        /*
        pos.x = pos.x / scale_.x;
        pos.y = pos.y / scale_.y;
        pos.z = pos.z / scale_.z;

        // world space -> uvw space
        float u = pos.x  / (float)Nx_;
        float v = pos.y  / (float)Ny_;
        float w = pos.z  / (float)Nz_;
        */

        // texture wraping
        // repeat
        u.x = u.x - floor(u.x);
        u.y = u.y - floor(u.y);
        u.z = u.z - floor(u.z);

        return u;
    }

    inline float sample_transmittance(const Float3& ws_pos) const {
        const Float3 uv_pos = uvw(ws_pos);
        const Float3 is_pos = times(uv_pos, Float3(Nx_, Ny_, Nz_));
        
        float ix = is_pos.x;
        float iy = is_pos.y;
        float iz = is_pos.z;
        /*
    
        int ix1 = (int)(ix + 0.5f) % Nx_;
        int ix0 = (ix1 - 1); if (ix0 < 0) ix0 += Nx_;
        float s0 = ix - (ix0 + 0.5f); if (s0 < 0) s0 += Nx_;
        float s1 = 1 - s0;

        int iy1 = (int)(iy + 0.5f) % Ny_;
        int iy0 = (iy1 - 1); if (iy0 < 0) iy0 += Ny_;
        float t0 = iy - (iy0 + 0.5f); if (t0 < 0) t0 += Ny_;
        float t1 = 1 - t0;

        int iz1 = (int)(iz + 0.5f) % Nz_;
        int iz0 = (iz1 - 1); if (iz0 < 0) iz0 += Nz_;
        float u0 = iz - (iz0 + 0.5f); if (u0 < 0) u0 += Nz_;
        float u1 = 1 - u0;
        
        const float density = 
            s0 * t0 * u0 * density_[I(ix0, iy0, iz0)] + 
            s0 * t0 * u1 * density_[I(ix0, iy0, iz1)] + 
            s0 * t1 * u0 * density_[I(ix0, iy1, iz0)] + 
            s0 * t1 * u1 * density_[I(ix0, iy1, iz1)] + 
            s1 * t0 * u0 * density_[I(ix1, iy0, iz0)] + 
            s1 * t0 * u1 * density_[I(ix1, iy0, iz1)] + 
            s1 * t1 * u0 * density_[I(ix1, iy1, iz0)] + 
            s1 * t1 * u1 * density_[I(ix1, iy1, iz1)];
        */
        const float density =  density_[I(ix, iy, iz)];

        return density / albedo_;
    }

    rt::BBox bbox_;
public:
    SimpleVolume(std::string filename, float density_scale, float albedo, Float3 &org, Float3 &scale) {
        FileManager fm;

        fm.load(filename.c_str());

        std::vector<unsigned char>& buf = fm.buffer();
        unsigned int *u_buf = (unsigned int*)&buf[0];


        Nx_ = u_buf[1];
        Ny_ = u_buf[2];
        Nz_ = u_buf[3];

        printf("Nx, Ny, Nz: %d %d %d\n", Nx_, Ny_, Nz_);

        density_.resize(Nx_ * Ny_ * Nz_);

        // float版
        float *f_buf = (float*)&u_buf[4];
        max_density_ = -1;
        for (int i = 0; i < Nx_ * Ny_ * Nz_; ++i) {
            density_[i] = density_scale * f_buf[i];

            if (max_density_ < density_[i])
                max_density_ = density_[i];
        }

        albedo_ = albedo;
        max_scattering_ = max_density_;
        max_transmittance_ = max_scattering_ / albedo_;

        org_ = org;
        scale_ = scale;
        bbox_ = rt::BBox(org_, org_ + times(scale_, Float3(Nx_, Ny_, Nz_)));

        factor_ = Float3(1 / scale_.x / Nx_, 1 / scale_.y / Ny_, 1 / scale_.z / Nz_);

        printf("Max Scattering: %f\n", max_scattering_);
    }

    bool inside(Float3& pt) const {
        return bbox_.inside(pt);
    }
    bool check_intersect(const rt::Ray& ray, float* hitt0, float* hitt1) const {
        return bbox_.check_intersect(ray, hitt0, hitt1);
    }
    
    float next_event(const rt::Ray& ray, Random& random) const {
        float u = 0;
        Float3 next_pos;

//		float d_max = hitt0 > 0 ? hitt0 : 1e32f;
        float d_max = 0;
        
        float hitt0 = -1, hitt1 = -1;

        if (bbox_.check_intersect(ray, &hitt0, &hitt1)) {
            d_max = std::max(hitt0, hitt1);
        }

        for (;;) {
            u += -log(random.nexto01()) / max_transmittance_;
            next_pos = ray.org + u * ray.dir;
            if (u >= d_max) {
                return -1.0f;
            }
            if (sample_transmittance(next_pos) / max_transmittance_ >= random.next01())
                break;
        }

        return u;
    }

};


Float3 radiance(Random& random, const rt::TriangleMesh& mesh, const SimpleVolume& volume, const rt::Ray& initial_ray) {

    const int kMaxScattering = 256;
    const float kAlbedo = 0.999f;

    rt::Ray current_ray = initial_ray;
    
    Color througput(1, 1, 1);
    Color L(0, 0, 0);
    for (int scattering = 0; scattering < kMaxScattering; ++scattering) {
        // VolumeのBBoxにヒットするか否か
        float hitt0 = -1, hitt1 = -1;
        bool is_volume_hit = false;
        if (volume.check_intersect(current_ray, &hitt0, &hitt1)) {
            is_volume_hit = true;
        }
        
        rt::Hitpoint current_hp;
        bool is_object_hit = false;
        if (!volume.inside(current_ray.org)) {
            if (is_volume_hit) {
                current_hp.distance = hitt0;
            }

            if (mesh.intersect_without_initialize(current_ray, &current_hp)) {
                is_object_hit = true;
                is_volume_hit = false;
            }
        }
        /*
        float distance_to_object = 1e+6f;
        rt::Hitpoint current_hp;
        bool is_object_hit = false;
        if (!volume.inside(current_ray.org) && mesh.intersect(current_ray, &current_hp)) {
            is_object_hit = true;
            distance_to_object = current_hp.distance;
        }

        // VolumeのBBoxにヒットするか否か
        float hitt0 = -1, hitt1 = -1;
        bool is_volume_hit = false;
        if (volume.check_intersect(current_ray, &hitt0, &hitt1)) {
            if (!is_object_hit || distance_to_object > hitt0) {
                is_volume_hit = true;
                is_object_hit = false;
            }
        }
        */

        if (!is_object_hit && !is_volume_hit) {
            L += Float3(0.01, 0.01, 0.01);
            break;
        }

        if (is_object_hit) {
            // 次の方向をBRDFインポータンスサンプリング
            const rt::TriangleElement t = mesh.getTriangle(current_hp.triangle_index);

            // 法線
            const Float3 e1 = normalize(*t.v[1] - *t.v[0]);
            const Float3 e2 = normalize(*t.v[2] - *t.v[0]);
            Float3 normal = normalize(cross(e1, e2));
            const Float3 hp_position = current_ray.org + current_hp.distance * current_ray.dir;

            const rt::Material* material = t.material;

            if (material == NULL)
                continue;

            rt::BRDF* current_brdf = NULL;
            rt::PhongBRDF phong(material->specular, material->specular_coefficient);
            rt::LambertianBRDF lambertian(material->diffuse);
                            
            if (random.next01() < material->metalic) {
                // スペキュラ
                current_brdf = &phong;
            } else {
                // ディフューズ
                current_brdf = &lambertian;
            }
            //current_brdf = &lambertian;

            Color color = current_brdf->reflectance();
            Float3 dir = current_brdf->sample(random, current_ray.dir, normal, NULL);
                            
            current_ray = rt::Ray(hp_position +  1e-2f * dir, dir);
            througput = times(color, througput);
        } else if (is_volume_hit) {
            const Float3 new_org0 = current_ray.org + (hitt0 + 1e-3f) * current_ray.dir;
            const Float3 new_org1 = current_ray.org + (hitt1 + 1e-3f) * current_ray.dir;
            if (!volume.inside(current_ray.org)) {
                // Volumeの外だったら中からスタート
                current_ray.org = new_org0;
            }

            const float tu = volume.next_event(current_ray, random);
            if (tu < 0) {
                // ボリューム外
                current_ray.org = new_org1;
                continue;
            }

            const Float3 next_pos = current_ray.org + tu * current_ray.dir;
            const bool cond2 = random.next01() < kAlbedo;
            if (!cond2) {
                througput = Float3(0, 0, 0);
                break;
            }
                    
            
            const int kLSample = 1;
            for (int i = 0; i < kLSample; ++i) {
                const Float3 next_dir = (Float3(0, 55 + 64, 0) / 50 - next_pos);
                const float phase = 1.0f / (4.0f * kPI);
                        
                        
                float u = volume.next_event(rt::Ray(next_pos, normalize(next_dir)), random);
                if (u < 0 || u >= next_dir.length()) {
                    L += phase * Float3(50, 2, 2) / kLSample;
                }
            }
            for (int i = 0; i < kLSample; ++i) {
                const Float3 next_dir = (Float3(15, -25 + 64, 8) / 50 - next_pos);
                const float phase = 1.0f / (4.0f * kPI);
                        
                float u = volume.next_event(rt::Ray(next_pos, normalize(next_dir)), random);
                if (u < 0 || u >= next_dir.length()) {
                    L += phase * Float3(2, 25, 2) / kLSample;
                }
                        
            }
            for (int i = 0; i < kLSample; ++i) {
                const Float3 next_dir = (Float3(-35, -5 + 64, -5) / 50 - next_pos);
                const float phase = 1.0f / (4.0f * kPI);
                        
                float u = volume.next_event(rt::Ray(next_pos, normalize(next_dir)), random);
                if (u < 0 || u >= next_dir.length()) {
                    L += phase * Float3(2, 2, 50) / kLSample;
                }
                        
            }


            // 位相関数でインポータンスサンプリング
            // 今回はとりあえず等方散乱としておく
            const Float3 next_dir = Sampling::uniformSphereSurface(random);
            const float phase = 1.0f / (4.0f * kPI);
            const float pdf = 1.0f / (4.0f * kPI);

            // throughtput = (phase / pdf) * throughtput;	
            current_ray = rt::Ray(next_pos, next_dir);
        }
    }

    return times(L, througput);
}

float to_bmp_value(const float value, const float display_gamma) {
    float f = (pow(value, 1.0f / display_gamma) * 255.0f + 0.5f);
    if (f >= 255)
        return 255;
    return f;
}


// imported from http://www.cs.rit.edu/~ncs/color/t_convert.html

// r,g,b values are from 0 to 1
// h = [0,360], s = [0,1], v = [0,1]
//		if s == 0, then h = -1 (undefined)
inline float MIN(float x, float y, float z) {
    return std::min(x, std::min(y, z));
}
inline float MAX(float x, float y, float z) {
    return std::max(x, std::max(y, z));
}

void RGBtoHSV( float r, float g, float b, float *h, float *s, float *v )
{
    float min, max, delta;
    min = MIN( r, g, b );
    max = MAX( r, g, b );
    *v = max;   // v
    delta = max - min;
    if( max != 0 )
        *s = delta / max;   // s
    else {
        // r = g = b = 0    // s = 0, v is undefined
        *s = 0;
        *h = -1;
        return;
    }
    if (delta == 0) {
        *h = 0;
        return;
    }

    if( r == max )
        *h = ( g - b ) / delta; // between yellow & magenta
    else if( g == max )
        *h = 2 + ( b - r ) / delta; // between cyan & yellow
    else
        *h = 4 + ( r - g ) / delta; // between magenta & cyan
    *h *= 60;   // degrees
    if( *h < 0 )
        *h += 360;
}
void HSVtoRGB( float *r, float *g, float *b, float h, float s, float v )
{
    int i;
    float f, p, q, t;
    if( s == 0 ) {
        // achromatic (grey)
        *r = *g = *b = v;
        return;
    }
    h /= 60;    // sector 0 to 5
    i = floor( h );
    f = h - i;  // factorial part of h
    p = v * ( 1 - s );
    q = v * ( 1 - s * f );
    t = v * ( 1 - s * ( 1 - f ) );
    switch( i ) {
        case 0:
            *r = v;
            *g = t;
            *b = p;
            break;
        case 1:
            *r = q;
            *g = v;
            *b = p;
            break;
        case 2:
            *r = p;
            *g = v;
            *b = t;
            break;
        case 3:
            *r = p;
            *g = q;
            *b = v;
            break;
        case 4:
            *r = t;
            *g = p;
            *b = v;
            break;
        default:    // case 5:
            *r = v;
            *g = p;
            *b = q;
            break;
    }
}

void rgb2YCoCg(float r, float g, float b, float *Y, float *Co, float *Cg) {
    *Y = 1/4.0f * r + 1/2.0f * g + 1/4.0f * b;
    *Co= 1/2.0f * r + 0/1.0f * g - 1/2.0f * b;
    *Cg=-1/4.0f * r + 1/2.0f * g - 1/4.0f * b;
}
void YCoCg2rgb(float Y, float Co, float Cg, float *r, float *g, float *b) {
    *r = Y + Co - Cg;
    *g = Y + 0  + Cg;
    *b = Y - Co - Cg;
}

inline void set_vector(std::vector<float>& arr, const int x, const int y, const int ch, const int width, const int height, float *vector, int learn_radius) {
    const int size = learn_radius * 2 + 1;

    for (int oy = -learn_radius; oy <= learn_radius; ++oy) {
        for (int ox = -learn_radius; ox <= learn_radius; ++ox) {
            const int nx = ox + x;
            const int ny = oy + y;
            if (0 <= nx && nx < width && 0 <= ny && ny < height) {
                vector[(oy + learn_radius) * size + (ox + learn_radius)] = arr[(ny * width + nx) * 3 + ch];
            }
        }
    }
}
inline float length(float *v0, float *v1, int size) {
    float sum = 0;
    for (int i = 0; i < size; ++i) {
        const float a = v0[i] - v1[i];
        sum += a * a;
    }
    return sum;
}

void nlm(std::vector<float>& arr0, std::vector<float>& arr1, const int width, const int height) {
    const int learn_radius = 3;
    const int compare_raidus = 6;
    const float sigma = 0.3f;

    // 変換
    for (int iy = 0; iy < height; ++iy) {
        for (int ix = 0; ix < width; ++ix) {
            const int idx = iy * width + ix;
            const float r = arr0[idx * 3 + 0];
            const float g = arr0[idx * 3 + 1];
            const float b = arr0[idx * 3 + 2];
            rgb2YCoCg(r, g, b, &arr0[idx * 3 + 0], &arr0[idx * 3 + 1], &arr0[idx * 3 + 2]);
        }
    }

    // NLM
    const int size = learn_radius * 2 + 1;
    for (int iy = 0; iy < height; ++iy) {
        std::cout << "Y: " << iy << "    \r";
        #pragma omp parallel for schedule(dynamic, 1)
        for (int ix = 0; ix < width; ++ix) {

            for (int ch = 0; ch < 3; ++ch) {
                float vector0[size * size] = {0};
                set_vector(arr0, ix, iy, ch, width, height, vector0, learn_radius);
                
                const int compare_size = compare_raidus * 2 + 1;
                float weight_map[compare_size * compare_size] = {0};
                float value_map[compare_size * compare_size] = {0};

                // 探索
                for (int oy = -compare_raidus; oy <= compare_raidus; ++oy) {
                    for (int ox = -compare_raidus; ox <= compare_raidus; ++ox) {
                        const int nx = ox + ix;
                        const int ny = oy + iy;
                        const int compare_idx = (oy + compare_raidus) * compare_size + (ox + compare_raidus);
                        if (0 <= nx && nx < width && 0 <= ny && ny < height) {
                            float vector1[size * size] = {0};
                            set_vector(arr0, nx, ny, ch, width, height, vector1, learn_radius);

                            // 重み計算
                            value_map[compare_idx] = arr0[(ny * width + nx) * 3 + ch];
                            weight_map[compare_idx] = length(vector0, vector1, size * size);
                        } else {
                            weight_map[compare_idx] = -1;
                        }
                    }
                }

                // 結果計算
                float sum = 0;
                float total_weight = 0;
                for (int cy = 0; cy < compare_size; ++cy) {
                    for (int cx = 0; cx < compare_size; ++ cx) {
                        const int compare_idx = cy * compare_size + cx;
                        if (weight_map[compare_idx] < 0)
                            continue;
                        const float weight = exp(-weight_map[compare_idx] / (sigma * sigma));
                        sum += value_map[compare_idx] * weight;
                        total_weight += weight;
                    }
                }
                if (total_weight > 0)
                    sum /= total_weight;
                arr1[(iy * width + ix) * 3 + ch] = sum;
            }
            const int idx = iy * width + ix;
            const float Y = arr1[idx * 3 + 0];
            const float Co= arr1[idx * 3 + 1];
            const float Cg= arr1[idx * 3 + 2];
            YCoCg2rgb(Y, Co, Cg, &arr1[idx * 3 + 0], &arr1[idx * 3 + 1], &arr1[idx * 3 + 2]);
        }
    }
}

/*
int main2() {
    using namespace hstd;
    Timer timer;
    timer.begin();
    Image image;
    HDROperator::load("I:/Projects/akari3/akari3/image_15.hdr", &image);

    const int width = image.width();
    const int height = image.height();
    
    std::vector<float> tmp_bmp_arr0(width * height * 3);
    std::vector<float> tmp_bmp_arr1(width * height * 3);
    for (int iy = 0; iy < height; ++iy) {
        for (int ix = 0; ix < width; ++ix) {
            const float display_gamma = 1.0f;

            float r = pow(image.at(ix, iy).x, 1 / display_gamma);
            float g = pow(image.at(ix, iy).y, 1 / display_gamma);
            float b = pow(image.at(ix, iy).z, 1 / display_gamma);

            tmp_bmp_arr0[(iy * width + ix) * 3 + 0] = r >= 1 ? 1 : r;
            tmp_bmp_arr0[(iy * width + ix) * 3 + 1] = g >= 1 ? 1 : g;
            tmp_bmp_arr0[(iy * width + ix) * 3 + 2] = b >= 1 ? 1 : b;
        }
    }

    nlm(tmp_bmp_arr0, tmp_bmp_arr1, width, height);

    
    for (int iy = 0; iy < height; ++iy) {
        for (int ix = 0; ix < width; ++ix) {
            const int idx = iy * width + ix;
            image.at(ix, iy) = Float3(tmp_bmp_arr1[idx * 3 + 0],  tmp_bmp_arr1[idx * 3 + 1], tmp_bmp_arr1[idx * 3 + 2]);  
        }
    }

    HDROperator::save("filtered.hdr", &image);
    std::cout << timer.end() / 1000.0f << " sec" << std::endl;
}
*/

int main(int argc, char *argv[]) {
    using namespace hstd;
    Timer timer;
    timer.begin();
    float accum_time = 0;
    float prev_now = timer.end();
    int image_iteration = 0;

    rt::TriangleMesh mesh;
    if (!rt::OBJOperator::load("asset/scene.obj", &mesh, Float3(1, 1, 1))) {
        return 0;
    }
   
    const int width = 1920 / 1;
    const int height = 1080 / 1;

    Image image(width, height);
    rt::PerspectiveCamera camera(Float3(-0.2, 1.3, 2.9), Float3(0.2, 0.9, 0), normalize(Float3(0.05, 1, 0)), 1, (float)height / width, 1.5f, 10000.f);
    //rt::PerspectiveCamera camera(Float3(1.2, 1.5, 1.5), Float3(0, 0.9, 0), Float3(0, 1, 0), 1, (float)height / width, 0.8f, 10000.f);
    
    const int kSample = 50 * 10;
    const int kMaxIteration = 1000000;

    const Color kSunPower = Color(30, 32, 33);
    const float kAlbedo = 0.999f;
    
    SimpleVolume ovolume("volume.float", 9.0f, kAlbedo,
        Float3(-128 * 0.02 / 2, 0.01, -128 * 0.02 / 2), Float3(0.02, 0.02, 0.02));
    
    std::vector<unsigned char> bmp_arr(width * height * 3);
    std::vector<float> tmp_bmp_arr0(width * height * 3);
    std::vector<float> tmp_bmp_arr1(width * height * 3);
    /*
    const float dist = (Float3(0, 0.9, 0) - Float3(2.6, 2.3, 5.1)).length();
    std::cout << dist << std::endl;
    */
    /*
    const float dist = (Float3(0, 0.9, 0) - Float3(0, 1.3, 3.1)).length();
    std::cout << dist << std::endl;
    */
    
    for (int iteration = 1; iteration < kMaxIteration; ++iteration) {
        std::cout << "iteration : " << iteration << "        \r";
        for (int iy = 0; iy < height; ++iy) {
            //std::cout << "y: " << iy << std::endl;
            #pragma omp parallel for schedule(dynamic, 1) num_threads(8)
            for (int ix = 0; ix < width; ++ix) {
                Random random((unsigned long long)(iy * width + ix) * kMaxIteration + iteration);

                //for (int ss = 0; ss < kSample; ++ss) {
                    const float u = (ix + random.next01()) / width;
                    const float v = (iy + random.next01()) / height;
                    rt::Ray ray = camera.get_ray(u, v);

                    // レンズ
                    float lensU, lensV;
                    const float r = random.next01();
                    const float lens_radius = 0.05f;
                    const float focal_distance = 5.89322; // 3.1257f;

                    Sampling::uniformCircle(random, &lensU, &lensV);
                    lensU *= lens_radius;
                    lensV *= lens_radius;
                    const float ft = fabs(focal_distance / dot(ray.dir, camera.camera_direction()));

                    const Float3 Pfocus = ray.org + ray.dir * ft;
                    ray.org = ray.org + lensU * camera.screen_side() +  lensV * camera.screen_up();
                    ray.dir = normalize(Pfocus - ray.org);

                    const Float3 L = radiance(random, mesh, ovolume, ray);

                    image.at(ix, iy) += L;
                //}
            }
        }

        // 出力チェック
        const float now = timer.end();

        float T = 0;
        if (image_iteration <= 28)
            T = 5.0f * 1000.0f;
        else
            T = 0.0f * 1000.0f; // nlm: 30sec

        if (accum_time >= T) {
            accum_time = 0;
            ++image_iteration;

            const float kScale = 1.3f;
        
            Image tmpimage(width, height);
            for (int iy = 0; iy < height; ++iy) {
                for (int ix = 0; ix < width; ++ix) {
                    tmpimage.at(ix, iy) = kScale * image.at(ix, iy) / iteration;
                }
            }

            // nlm
            const float display_gamma = 2.2f;

            // 最後だけ
            if (image_iteration == 30) {
                for (int iy = 0; iy < height; ++iy) {
                    for (int ix = 0; ix < width; ++ix) {
                    
                        float r = pow(tmpimage.at(ix, iy).x, 1 / display_gamma);
                        float g = pow(tmpimage.at(ix, iy).y, 1 / display_gamma);
                        float b = pow(tmpimage.at(ix, iy).z, 1 / display_gamma);

                        tmp_bmp_arr0[(iy * width + ix) * 3 + 0] = b >= 1 ? 1 : b;
                        tmp_bmp_arr0[(iy * width + ix) * 3 + 1] = g >= 1 ? 1 : g;
                        tmp_bmp_arr0[(iy * width + ix) * 3 + 2] = r >= 1 ? 1 : r;
                    }
                }

                nlm(tmp_bmp_arr0, tmp_bmp_arr1, width, height);
                std::cout << std::endl;
    
                for (int iy = 0; iy < height; ++iy) {
                    for (int ix = 0; ix < width; ++ix) {
                        const int idx = iy * width + ix;
                        tmpimage.at(ix, iy) = Float3(tmp_bmp_arr1[idx * 3 + 0],  tmp_bmp_arr1[idx * 3 + 1], tmp_bmp_arr1[idx * 3 + 2]);  
                    }
                }    
                // bmpに出力
                for (int iy = 0; iy < height; ++iy) {
                    for (int ix = 0; ix < width; ++ix) {
                        float r = tmpimage.at(ix, iy).x * 255.0f;
                        float g = tmpimage.at(ix, iy).y * 255.0f;
                        float b = tmpimage.at(ix, iy).z * 255.0f;

                        bmp_arr[(iy * width + ix) * 3 + 0] = r >= 255 ? 255 : r;
                        bmp_arr[(iy * width + ix) * 3 + 1] = g >= 255 ? 255 : g;
                        bmp_arr[(iy * width + ix) * 3 + 2] = b >= 255 ? 255 : b;
                    }
                }
            } else {
                // bmpに出力
                for (int iy = 0; iy < height; ++iy) {
                    for (int ix = 0; ix < width; ++ix) {
                        float r = to_bmp_value(tmpimage.at(ix, iy).x, display_gamma);
                        float g = to_bmp_value(tmpimage.at(ix, iy).y, display_gamma);
                        float b = to_bmp_value(tmpimage.at(ix, iy).z, display_gamma);

                        bmp_arr[(iy * width + ix) * 3 + 0] = b >= 255 ? 255 : b;
                        bmp_arr[(iy * width + ix) * 3 + 1] = g >= 255 ? 255 : g;
                        bmp_arr[(iy * width + ix) * 3 + 2] = r >= 255 ? 255 : r;
                    }
                }
            }
            


            char buf[1024];
            sprintf(buf, "%03d.bmp", image_iteration);
            exportToBmp(buf, &bmp_arr[0], width, height);
            /*
            sprintf(buf, "image_%02d.hdr", image_iteration);
            HDROperator::save(buf, &tmpimage, true);
            */
            std::cout << std::endl << "Save: " << buf << std::endl << std::endl;

            if (image_iteration == 30) {
                std::cout << timer.end() / 1000.0f << " sec" << std::endl;
                return 0;
            }
        } else {
            accum_time += now - prev_now;
        }
        prev_now = now;
    }

    
    // HDROperator::save("result.hdr", &image);
}