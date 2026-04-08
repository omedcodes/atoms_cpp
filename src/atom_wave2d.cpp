#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <type_traits>
#include <vector>
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif
using namespace glm;
using namespace std;

float oribtDistance = 200.0f;

vec2 mouseWorld(0.0f);
struct Engine {
    GLFWwindow* window;
    int WIDTH = 1200;
    int HEIGHT = 800;

    Engine(){
        if (!glfwInit()){
            cerr << "failed to initialize glfw, HOW-";
            exit(EXIT_FAILURE);
        }

        window = glfwCreateWindow(WIDTH, HEIGHT, "2D Simulation of Atoms", nullptr, nullptr);
        if(!window){
            cerr << "failed to create window, WDYM?!";
            glfwTerminate();
            exit(EXIT_FAILURE);
        }

        glfwMakeContextCurrent(window);
        int fbWidth, fbHeight;
        glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
        glViewport(0, 0, fbWidth, fbHeight);
    }
    void run(){
        glClear(GL_COLOR_BUFFER_BIT);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        // find origin at center
        double halfWidth = WIDTH / 2.0;
        double halfHeight = HEIGHT / 2.0;

        glOrtho(-halfWidth, halfWidth, -halfHeight, halfHeight, -1.0, 1.0);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    }
};
Engine engine;
