#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#include <cmath>
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

int N = 10000;
float LightingScaler = 800;
float n = 3; float l = 1; float m = 1;

struct Particle {
    vec3 pos;
    vec3 vel = vec3(0.0f);
    vec4 color;
    Particle(vec3 p, vec4 c = vec4(0.0f, 0.5f, 1.0f, 1.0f)) : pos(p), color(c){}
};
vector<Particle> particles;

random_device rd; mt19937 gen(rd()); uniform_real_distribution<float> dis(0.0f, 1.0f);

struct CDFCache {
    vector<double> rCDF, thetaCDF;
    int cached_n = -1, cached_l = -1, cached_m = -1;
    bool rValid = false, thetaValid = false;
} cdfCache;

double sampleR(int n, int l, mt19937& gen) {
    const int NCDF = 4096;
    const double rMax = 10.0 * n * n * a0;

    if (!cdfCache.rValid || cdfCache.cached_n != n || cdfCache.cached_l != l) {
        cdfCache.rCDF.resize(NCDF);
        double sum = 0.0;
        for (int i = 0; i < NCDF; ++i) {
            double r = i * (rMax / (NCDF - 1));
            double rho = 2.0 * r / (n * a0);
            int k = n - l - 1;
            int alpha = 2 * l + 1;
            double L = 1.0, Lm1 = 1.0 + alpha - rho;
            if (k == 1) L = Lm1;
            else if (k > 1) {
                double Lm2 = 1.0;
                for (int j = 2; j <= k; ++j) {
                    L = ((2*j - 1 + alpha - rho) * Lm1 - (j - 1 + alpha) * Lm2) / j;
                    Lm2 = Lm1; Lm1 = L;
                }
            }
            double norm = pow(2.0 / (n * a0), 3) * tgamma(n - l) / (2.0 * n * tgamma(n + l + 1));
            double R = sqrt(norm) * exp(-rho / 2.0) * pow(rho, l) * L;
            sum += r * r * R * R;
            cdfCache.rCDF[i] = sum;
        }
        for (double& v : cdfCache.rCDF) v /= sum;
        cdfCache.cached_n = n; cdfCache.cached_l = l;
        cdfCache.rValid = true;
    }

    uniform_real_distribution<double> d(0.0, 1.0);
    int idx = (int)(lower_bound(cdfCache.rCDF.begin(), cdfCache.rCDF.end(), d(gen)) - cdfCache.rCDF.begin());
    return idx * (10.0 * n * n * a0 / 4095.0);
}

double sampleTheta(int l, int m, mt19937& gen) {
    const int NCDF = 2048;

    if (!cdfCache.thetaValid || cdfCache.cached_l != l || cdfCache.cached_m != (int)m) {
        cdfCache.thetaCDF.resize(NCDF);
        double sum = 0.0;
        for (int i = 0; i < NCDF; ++i) {
            double theta = i * (M_PI / (NCDF - 1));
            double x = cos(theta);
            double Pmm = 1.0;
            if (m > 0) {
                double somx2 = sqrt((1.0 - x) * (1.0 + x));
                double fact = 1.0;
                for (int j = 1; j <= m; ++j) { Pmm *= -fact * somx2; fact += 2.0; }
            }
            double Plm;
            if (l == m) { Plm = Pmm; }
            else {
                double Pm1m = x * (2 * m + 1) * Pmm;
                if (l == m + 1) { Plm = Pm1m; }
                else {
                    double Pll;
                    for (int ll = m + 2; ll <= l; ++ll) {
                        Pll = ((2 * ll - 1) * x * Pm1m - (ll + m - 1) * Pmm) / (ll - m);
                        Pmm = Pm1m; Pm1m = Pll;
                    }
                    Plm = Pm1m;
                }
            }
            sum += sin(theta) * Plm * Plm;
            cdfCache.thetaCDF[i] = sum;
        }
        for (double& v : cdfCache.thetaCDF) v /= sum;
        cdfCache.cached_l = l; cdfCache.cached_m = (int)m;
        cdfCache.thetaValid = true;
    }

    uniform_real_distribution<double> d(0.0, 1.0);
    int idx = (int)(lower_bound(cdfCache.thetaCDF.begin(), cdfCache.thetaCDF.end(), d(gen)) - cdfCache.thetaCDF.begin());
    return idx * (M_PI / 2047.0);
}

float samplePhi() {
    return 2.0f * M_PI * dis(gen);
}

vec3 calculateProbabilityFlow(Particle& p, int n, int l, int m) {
    double r = length(p.pos); if (r < 1e-6) return vec3(0.0f);
    double theta = acos(p.pos.y / r);
    double phi = atan2(p.pos.z, p.pos.x);
    double sinTheta = sin(theta); if (abs(sinTheta) < 1e-4) sinTheta = 1e-4;
    double v_mag = hbar * m / (m_e * r * sinTheta);
    return vec3((float)(-v_mag * sin(phi)), 0.0f, (float)(v_mag * cos(phi)));
}

vec4 heatmap_fire(float value) {
    value = std::max(0.0f, std::min(1.0f, value));
    const int num_stops = 6;
    vec4 colors[num_stops] = {
        {0.0f, 0.0f, 0.0f, 1.0f},
        {0.3f, 0.0f, 0.6f, 1.0f},
        {0.8f, 0.0f, 0.0f, 1.0f},
        {1.0f, 0.5f, 0.0f, 1.0f},
        {1.0f, 1.0f, 0.0f, 1.0f},
        {1.0f, 1.0f, 1.0f, 1.0f}
    };
    float sv = value * (num_stops - 1);
    int i = std::min((int)sv, num_stops - 2);
    float t = sv - i;
    return vec4(
        colors[i].r + t * (colors[i+1].r - colors[i].r),
        colors[i].g + t * (colors[i+1].g - colors[i].g),
        colors[i].b + t * (colors[i+1].b - colors[i].b),
        1.0f
    );
}

vec4 inferno(double r, double theta, int n, int l, int m) {
    double rho = 2.0 * r / (n * a0);
    int k = n - l - 1;
    int alpha = 2 * l + 1;
    double L = 1.0;
    if (k == 1) { L = 1.0 + alpha - rho; }
    else if (k > 1) {
        double Lm2 = 1.0, Lm1 = 1.0 + alpha - rho;
        for (int j = 2; j <= k; ++j) {
            L = ((2*j - 1 + alpha - rho) * Lm1 - (j - 1 + alpha) * Lm2) / j;
            Lm2 = Lm1; Lm1 = L;
        }
    }
    double norm = pow(2.0 / (n * a0), 3) * tgamma(n - l) / (2.0 * n * tgamma(n + l + 1));
    double R = sqrt(norm) * exp(-rho / 2.0) * pow(rho, l) * L;
    double x = cos(theta);
    double Pmm = 1.0;
    if (m > 0) {
        double somx2 = sqrt((1.0 - x) * (1.0 + x));
        double fact = 1.0;
        for (int j = 1; j <= m; ++j) { Pmm *= -fact * somx2; fact += 2.0; }
    }
    double Plm;
    if (l == m) { Plm = Pmm; }
    else {
        double Pm1m = x * (2*m + 1) * Pmm;
        if (l == m + 1) { Plm = Pm1m; }
        else {
            for (int ll = m + 2; ll <= l; ++ll) {
                double Pll = ((2*ll - 1) * x * Pm1m - (ll + m - 1) * Pmm) / (ll - m);
                Pmm = Pm1m; Pm1m = Pll;
            }
            Plm = Pm1m;
        }
    }
    return heatmap_fire((R * R) * (Plm * Plm) * LightingScaler);
}

struct Camera {
    vec3 target = vec3(0.0f);
    float radius = 50.0f;
    float azimuth = 0.0f;
    float elevation = M_PI / 2.0f;
    float orbitSpeed = 0.01f;
    double zoomSpeed = zmSpeed;
    bool dragging = false;
    double lastX = 0.0, lastY = 0.0;

    vec3 position() const {
        float ce = clamp(elevation, 0.01f, float(M_PI) - 0.01f);
        return vec3(radius * sin(ce) * cos(azimuth), radius * cos(ce), radius * sin(ce) * sin(azimuth));
    }
    void processMouseMove(double x, double y) {
        if (dragging) {
            azimuth += float(x - lastX) * orbitSpeed;
            elevation = clamp(elevation - float(y - lastY) * orbitSpeed, 0.01f, float(M_PI) - 0.01f);
        }
        lastX = x; lastY = y;
    }
    void processMouseButton(int button, int action, int mods, GLFWwindow* win) {
        if (button == GLFW_MOUSE_BUTTON_LEFT || button == GLFW_MOUSE_BUTTON_MIDDLE) {
            if (action == GLFW_PRESS) { dragging = true; glfwGetCursorPos(win, &lastX, &lastY); }
            else if (action == GLFW_RELEASE) dragging = false;
        }
    }
    void processScroll(double xoffset, double yoffset) {
        radius = std::max(1.0f, radius - (float)(yoffset * zoomSpeed));
    }
};
Camera camera;

vec3 sphericalToCartesian(float r, float theta, float phi){
    return vec3(r * sin(theta) * cos(phi), r * cos(theta), r * sin(theta) * sin(phi));
}

void generateParticles(int N) {
    cdfCache.rValid = false;
    cdfCache.thetaValid = false;
    particles.clear();
    particles.reserve(N);
    for (int i = 0; i < N; ++i) {
        float r     = sampleR(n, l, gen);
        float theta = sampleTheta(l, m, gen);
        float phi   = samplePhi();
        vec3 pos = sphericalToCartesian(r, theta, phi);
        float pr = length(pos);
        vec4 col = inferno(pr, acos(pos.y / std::max(pr, 1e-6f)), n, l, m);
        particles.emplace_back(pos, col);
    }
}

struct InstanceData {
    vec3 pos;
    float pad;
    vec4 color;
};

struct Engine {
    GLFWwindow* window;
    int WIDTH = 800;
    int HEIGHT = 600;

    GLuint instancedShaderProgram;
    GLuint sphereVAO, sphereVBO, sphereEBO;
    GLuint instanceVBO;
    int sphereIndexCount = 0;

    void buildSphereGeometry(int stacks, int slices) {
        vector<float> verts;
        vector<unsigned int> indices;
        for (int i = 0; i <= stacks; ++i) {
            float phi = M_PI * i / stacks;
            for (int j = 0; j <= slices; ++j) {
                float theta = 2.0f * M_PI * j / slices;
                verts.push_back(sin(phi) * cos(theta));
                verts.push_back(cos(phi));
                verts.push_back(sin(phi) * sin(theta));
            }
        }
        for (int i = 0; i < stacks; ++i) {
            for (int j = 0; j < slices; ++j) {
                unsigned int a = i * (slices + 1) + j;
                unsigned int b = a + slices + 1;
                indices.push_back(a); indices.push_back(b); indices.push_back(a + 1);
                indices.push_back(b); indices.push_back(b + 1); indices.push_back(a + 1);
            }
        }
        sphereIndexCount = (int)indices.size();

        glGenVertexArrays(1, &sphereVAO);
        glGenBuffers(1, &sphereVBO);
        glGenBuffers(1, &sphereEBO);
        glBindVertexArray(sphereVAO);
        glBindBuffer(GL_ARRAY_BUFFER, sphereVBO);
        glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(float), verts.data(), GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphereEBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

        glGenBuffers(1, &instanceVBO);
        glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
        glBufferData(GL_ARRAY_BUFFER, N * sizeof(InstanceData), nullptr, GL_DYNAMIC_DRAW);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (void*)offsetof(InstanceData, pos));
        glVertexAttribDivisor(1, 1);
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (void*)offsetof(InstanceData, color));
        glVertexAttribDivisor(2, 1);
        glBindVertexArray(0);
    }

    GLuint compileShader(GLenum type, const char* src) {
        GLuint s = glCreateShader(type);
        glShaderSource(s, 1, &src, NULL);
        glCompileShader(s);
        GLint ok; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
        if (!ok) { char log[512]; glGetShaderInfoLog(s, 512, NULL, log); cerr << log << "\n"; }
        return s;
    }

    GLuint buildInstancedProgram() {
        const char* vert = R"(
#version 330 core
layout(location = 0) in vec3 aLocalPos;
layout(location = 1) in vec3 aInstancePos;
layout(location = 2) in vec4 aColor;
uniform mat4 uVP;
uniform float uRadius;
uniform vec3 uLightPos;
out vec4 vColor;
void main() {
    vec3 worldPos = aLocalPos * uRadius + aInstancePos;
    vec3 lightDir = normalize(uLightPos - worldPos);
    float diff = max(dot(aLocalPos, lightDir), 0.0);
    vColor = vec4(aColor.rgb * (0.25 + 0.75 * diff), aColor.a);
    gl_Position = uVP * vec4(worldPos, 1.0);
}
)";
        const char* frag = R"(
#version 330 core
in vec4 vColor;
out vec4 FragColor;
void main() { FragColor = vColor; }
)";
        GLuint vs = compileShader(GL_VERTEX_SHADER, vert);
        GLuint fs = compileShader(GL_FRAGMENT_SHADER, frag);
        GLuint prog = glCreateProgram();
        glAttachShader(prog, vs); glAttachShader(prog, fs);
        glLinkProgram(prog);
        glDeleteShader(vs); glDeleteShader(fs);
        return prog;
    }

    Engine() {
        if (!glfwInit()) { cerr << "GLFW init failed\n"; exit(EXIT_FAILURE); }
        window = glfwCreateWindow(WIDTH, HEIGHT, "Hydrogen Atom", nullptr, nullptr);
        if (!window) { cerr << "Failed to create window\n"; glfwTerminate(); exit(EXIT_FAILURE); }
        glfwMakeContextCurrent(window);
        glfwSwapInterval(1);
        glViewport(0, 0, WIDTH, HEIGHT);
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glewExperimental = GL_TRUE;
        if (glewInit() != GLEW_OK) { cerr << "GLEW init failed\n"; glfwTerminate(); exit(EXIT_FAILURE); }
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        instancedShaderProgram = buildInstancedProgram();
        buildSphereGeometry(6, 8);
    }

    void onResize(int w, int h) {
        if (w == 0 || h == 0) return;
        WIDTH = w; HEIGHT = h;
        glViewport(0, 0, w, h);
    }

    void resizeInstanceBuffer(int newN) {
        glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
        glBufferData(GL_ARRAY_BUFFER, newN * sizeof(InstanceData), nullptr, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    vector<InstanceData> instanceStagingBuf;

    void render() {
        instanceStagingBuf.resize(particles.size());
        for (int i = 0; i < (int)particles.size(); ++i) {
            instanceStagingBuf[i].pos   = particles[i].pos;
            instanceStagingBuf[i].pad   = 0.0f;
            instanceStagingBuf[i].color = particles[i].color;
        }
        glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, instanceStagingBuf.size() * sizeof(InstanceData), instanceStagingBuf.data());
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glUseProgram(instancedShaderProgram);

        mat4 view = lookAt(camera.position(), camera.target, vec3(0,1,0));
        mat4 proj = perspective(radians(45.0f), (float)WIDTH / HEIGHT, 0.1f, 10000.0f);
        mat4 vp = proj * view;

        glUniformMatrix4fv(glGetUniformLocation(instancedShaderProgram, "uVP"),      1, GL_FALSE, value_ptr(vp));
        glUniform1f(glGetUniformLocation(instancedShaderProgram, "uRadius"),  electron_r);
        glUniform3fv(glGetUniformLocation(instancedShaderProgram, "uLightPos"),1, value_ptr(vec3(0.0f, 50.0f, 50.0f)));

        glBindVertexArray(sphereVAO);
        glDrawElementsInstanced(GL_TRIANGLES, sphereIndexCount, GL_UNSIGNED_INT, 0, (GLsizei)particles.size());
        glBindVertexArray(0);
    }

    void setupCallbacks() {
        glfwSetWindowUserPointer(window, this);

        glfwSetFramebufferSizeCallback(window, [](GLFWwindow* win, int w, int h) {
            ((Engine*)glfwGetWindowUserPointer(win))->onResize(w, h);
        });
        glfwSetMouseButtonCallback(window, [](GLFWwindow* win, int button, int action, int mods) {
            camera.processMouseButton(button, action, mods, win);
        });
        glfwSetCursorPosCallback(window, [](GLFWwindow* win, double x, double y) {
            camera.processMouseMove(x, y);
        });
        glfwSetScrollCallback(window, [](GLFWwindow* win, double xoffset, double yoffset) {
            camera.processScroll(xoffset, yoffset);
        });
        glfwSetKeyCallback(window, [](GLFWwindow* win, int key, int scancode, int action, int mods) {
            if (action != GLFW_PRESS) return;
            Engine* eng = (Engine*)glfwGetWindowUserPointer(win);

            if (key == GLFW_KEY_W) { n += 1; }
            else if (key == GLFW_KEY_S) { n -= 1; if (n < 1) n = 1; }
            else if (key == GLFW_KEY_E) { l += 1; }
            else if (key == GLFW_KEY_D) { l -= 1; if (l < 0) l = 0; }
            else if (key == GLFW_KEY_R) { m += 1; }
            else if (key == GLFW_KEY_F) { m -= 1; }
            else if (key == GLFW_KEY_T) { N = std::min(N * 2, 50000); eng->resizeInstanceBuffer(N); }
            else if (key == GLFW_KEY_ESCAPE) { glfwSetWindowShouldClose(win, GLFW_TRUE); return; }

            if (l > n - 1) l = n - 1;
            if (l < 0) l = 0;
            if (m > l) m = l;
            if (m < -l) m = -l;

            generateParticles(N);
            cout << "n=" << n << " l=" << l << " m=" << m << " N=" << N << "\n";
        });
    }
};

int main() {
    Engine engine;
    engine.setupCallbacks();
    generateParticles(N);

    float dt = 0.5f;

    double lastTime = glfwGetTime();
    int frameCount = 0;

    while (!glfwWindowShouldClose(engine.window)) {
        for (Particle& p : particles) {
            double r = length(p.pos);
            if (r > 1e-6) {
                double theta  = acos(p.pos.y / r);
                double new_phi = atan2(p.pos.z + calculateProbabilityFlow(p, n, l, m).z * dt,
                                       p.pos.x + calculateProbabilityFlow(p, n, l, m).x * dt);
                p.pos = sphericalToCartesian(r, theta, new_phi);
            }
        }

        engine.render();
        glfwSwapBuffers(engine.window);
        glfwPollEvents();

        ++frameCount;
        double now = glfwGetTime();
        if (now - lastTime >= 1.0) {
            cout << "FPS: " << frameCount << "\n";
            frameCount = 0;
            lastTime = now;
        }
    }

    glfwDestroyWindow(engine.window);
    glfwTerminate();
    return 0;
}