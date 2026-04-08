#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <thread>
#include <chrono>
#include <fstream>
#include <complex>
#include <random>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace glm;
using namespace std;

const float a0 = 1;
const float electron_r = 0.25f;
const double hbar = 1;
const double m_e = 1;
const double zmSpeed = 10.0;

int N = 100000;
float LightingScaler = 700;
float n = 3; float l = 1; float m = 1;

struct Particle {
    vec3 pos;
    vec3 vel = vec3(0.0f);
    vec4 color;
    Particle(vec3 p, vec4 c = vec4(0.0f, 0.5f, 1.0f, 1.0f)) : pos(p), color(c){}
};
vector<Particle> particles;

random_device rd; mt19937 gen(rd()); uniform_real_distribution<float> dis(0.0f, 1.0f);

double sampleR(int n, int l, mt19937& gen) {
    const int N = 4096;
    //const double a0 = 1.0;
    const double rMax = 10.0 * n * n * a0;

    static vector<double> cdf;
    static bool built = false;

    if (!built) {
        cdf.resize(N);
        double dr = rMax / (N - 1);
        double sum = 0.0;

        for (int i = 0; i < N; ++i) {
            double r = i * dr;
            double rho = 2.0 * r / (n * a0);

            int k = n - l - 1;
            int alpha = 2 * l + 1;

            double L = 1.0, Lm1 = 1.0 + alpha - rho;
            if (k == 1) L = Lm1;
            else if (k > 1) {
                double Lm2 = 1.0;
                for (int j = 2; j <= k; ++j) {
                    L = ((2*j - 1 + alpha - rho) * Lm1 -
                         (j - 1 + alpha) * Lm2) / j;
                    Lm2 = Lm1;
                    Lm1 = L;
                }
            }

            double norm = pow(2.0 / (n * a0), 3) * tgamma(n - l) / (2.0 * n * tgamma(n + l + 1));
            double R = sqrt(norm) * exp(-rho / 2.0) * pow(rho, l) * L;

            double pdf = r * r * R * R;
            sum += pdf;
            cdf[i] = sum;
        }

        for (double& v : cdf) v /= sum;
        built = true;
    }

    uniform_real_distribution<double> dis(0.0, 1.0);
    double u = dis(gen);

    int idx = lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin();
    return idx * (rMax / (N - 1));
}
double sampleTheta(int l, int m, mt19937& gen) {
    const int N = 2048;
    static vector<double> cdf;
    static bool built = false;

    if (!built) {
        cdf.resize(N);
        double dtheta = M_PI / (N - 1);
        double sum = 0.0;

        for (int i = 0; i < N; ++i) {
            double theta = i * dtheta;
            double x = cos(theta);

            double Pmm = 1.0;
            if (m > 0) {
                double somx2 = sqrt((1.0 - x) * (1.0 + x));
                double fact = 1.0;
                for (int j = 1; j <= m; ++j) {
                    Pmm *= -fact * somx2;
                    fact += 2.0;
                }
            }

            double Plm;
            if (l == m) {
                Plm = Pmm;
            } else {
                double Pm1m = x * (2 * m + 1) * Pmm;
                if (l == m + 1) {
                    Plm = Pm1m;
                } else {
                    double Pll;
                    for (int ll = m + 2; ll <= l; ++ll) {
                        Pll = ((2 * ll - 1) * x * Pm1m -
                               (ll + m - 1) * Pmm) / (ll - m);
                        Pmm = Pm1m;
                        Pm1m = Pll;
                    }
                    Plm = Pm1m;
                }
            }

            double pdf = sin(theta) * Plm * Plm;
            sum += pdf;
            cdf[i] = sum;
        }

        for (double& v : cdf) v /= sum;
        built = true;
    }

    uniform_real_distribution<double> dis(0.0, 1.0);
    double u = dis(gen);

    int idx = lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin();
    return idx * (M_PI / (N - 1));
}
float samplePhi(float n, float l, float m) {
    return 2.0f * M_PI * dis(gen);
}
vec3 calculateProbabilityFlow(Particle& p, int n, int l, int m) {
    double r = length(p.pos);   if (r < 1e-6) return vec3(0.0f);
    double theta = acos(p.pos.y / r); 
    double phi = atan2(p.pos.z, p.pos.x); 

    double sinTheta = sin(theta);  if (abs(sinTheta) < 1e-4) sinTheta = 1e-4;
    double v_mag = hbar * m / (m_e * r * sinTheta);

    double vx = -v_mag * sin(phi);
    double vy = 0.0; 
    double vz =  v_mag * cos(phi);

    return vec3((float)vx, (float)vy, (float)vz);
}

vec4 inferno2(double r, double theta, double phi, int n, int l, int m)
{
    double rho = 2.0 * r / (n * a0);

    int k = n - l - 1;
    int alpha = 2 * l + 1;

    double L = 1.0;
    if (k == 1) {
        L = 1.0 + alpha - rho;
    } else if (k > 1) {
        double Lm2 = 1.0;
        double Lm1 = 1.0 + alpha - rho;
        for (int j = 2; j <= k; ++j) {
            L = ((2*j - 1 + alpha - rho) * Lm1 -
                 (j - 1 + alpha) * Lm2) / j;
            Lm2 = Lm1;
            Lm1 = L;
        }
    }

    double norm = pow(2.0 / (n * a0), 3) * tgamma(n - l) / (2.0 * n * tgamma(n + l + 1));

    double R = sqrt(norm) * exp(-rho / 2.0) * pow(rho, l) * L;
    double radial = R * R;
    double x = cos(theta);

    double Pmm = 1.0;
    if (m > 0) {
        double somx2 = sqrt((1.0 - x) * (1.0 + x));
        double fact = 1.0;
        for (int j = 1; j <= m; ++j) {
            Pmm *= -fact * somx2;
            fact += 2.0;
        }
    }

    double Plm;
    if (l == m) {
        Plm = Pmm;
    } else {
        double Pm1m = x * (2*m + 1) * Pmm;
        if (l == m + 1) {
            Plm = Pm1m;
        } else {
            for (int ll = m + 2; ll <= l; ++ll) {
                double Pll = ((2*ll - 1) * x * Pm1m - (ll + m - 1) * Pmm) / (ll - m);
                Pmm = Pm1m;
                Pm1m = Pll;
            }
            Plm = Pm1m;
        }
    }

    double angular = Plm * Plm;
    double intensity = radial * angular;

    double t = log10(intensity + 1e-12) + 12.0;
    t /= 12.0;

    t = clamp(t, 0.0, 1.0);

    float rC = smoothstep(0.15f, 1.0f, static_cast<float>(t));
    float gC = smoothstep(0.45f, 1.0f, static_cast<float>(t));
    float bC = smoothstep(0.85f, 1.0f, static_cast<float>(t)) * 0.2f;

    return vec4(rC, gC * 0.8f, bC, 1.0f);
}

vec4 heatmap_fire(float value) {
    value = std::max(0.0f, std::min(1.0f, value));

    // Define color stops for the "Heat/Fire" pattern
    // Order: Black -> Dark Purple -> Red -> Orange -> Yellow -> White
    const int num_stops = 6;
    vec4 colors[num_stops] = {
        {0.0f, 0.0f, 0.0f, 1.0f}, // 0.0: Black
        {0.3f, 0.0f, 0.6f, 1.0f}, // 0.2: Dark Purple
        {0.8f, 0.0f, 0.0f, 1.0f}, // 0.4: Deep Red
        {1.0f, 0.5f, 0.0f, 1.0f}, // 0.6: Orange
        {1.0f, 1.0f, 0.0f, 1.0f}, // 0.8: Yellow
        {1.0f, 1.0f, 1.0f, 1.0f}  // 1.0: White
    };

    // Find which segment the value falls into
    float scaled_v = value * (num_stops - 1);
    int i = static_cast<int>(scaled_v);
    int next_i = std::min(i + 1, num_stops - 1);
    
    // Calculate how far we are between stop 'i' and 'next_i'
    float local_t = scaled_v - i;

    vec4 result;
    result.r = colors[i].r + local_t * (colors[next_i].r - colors[i].r);
    result.g = colors[i].g + local_t * (colors[next_i].g - colors[i].g);
    result.b = colors[i].b + local_t * (colors[next_i].b - colors[i].b);
    result.a = 1.0f; // Solid opacity
    // result = vec4(0.2, 0.9, 0.05, 1.0);

    return result;
}
vec4 inferno(double r, double theta, double phi, int n, int l, int m) {
    // radial part |R(r)|^2 
    double rho = 2.0 * r / (n * a0);

    int k = n - l - 1;
    int alpha = 2 * l + 1;

    double L = 1.0;
    if (k == 1) {
        L = 1.0 + alpha - rho;
    } else if (k > 1) {
        double Lm2 = 1.0;
        double Lm1 = 1.0 + alpha - rho;
        for (int j = 2; j <= k; ++j) {
            L = ((2*j - 1 + alpha - rho) * Lm1 -
                 (j - 1 + alpha) * Lm2) / j;
            Lm2 = Lm1;
            Lm1 = L;
        }
    }

    double norm = pow(2.0 / (n * a0), 3) * tgamma(n - l) / (2.0 * n * tgamma(n + l + 1));

    double R = sqrt(norm) * exp(-rho / 2.0) * pow(rho, l) * L;
    double radial = R * R;

    // angular part |P_l^m(cosθ)|^2 
    double x = cos(theta);

    double Pmm = 1.0;
    if (m > 0) {
        double somx2 = sqrt((1.0 - x) * (1.0 + x));
        double fact = 1.0;
        for (int j = 1; j <= m; ++j) {
            Pmm *= -fact * somx2;
            fact += 2.0;
        }
    }

    double Plm;
    if (l == m) {
        Plm = Pmm;
    } else {
        double Pm1m = x * (2*m + 1) * Pmm;
        if (l == m + 1) {
            Plm = Pm1m;
        } else {
            for (int ll = m + 2; ll <= l; ++ll) {
                double Pll = ((2*ll - 1) * x * Pm1m - (ll + m - 1) * Pmm) / (ll - m);
                Pmm = Pm1m;
                Pm1m = Pll;
            }
            Plm = Pm1m;
        }
    }

    double angular = Plm * Plm;
    double intensity = radial * angular;

    //cout << "intensity: " << intensity << endl;
    // return vec4(1.0f);
    return heatmap_fire(intensity * LightingScaler);
}
