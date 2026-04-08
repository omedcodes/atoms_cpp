# Atoms Cpp

This project contains high performance simulations of hydrogen wave functions and raytracing. The code is written in C++ and utilizes OpenGL for the rendering.

## Build

You need a compiler like Clang or G++ and the development files for GLFW and GLEW.

On Windows use `mingw32-make`
On Linux or Mac use `make`

The build system creates a bin directory where the executables are stored. The Windows build is configured to link everything statically so the exe files can be shared and run without additional setup or dll files.

## Usage

Once compiled you can run the programs directly from the bin folder. 

Atom Raytracer: Visualizes the probability density using raytracing.
Atom Wave2D: A two dimensional representation of the wave function flow.

Movement is usually handled via mouse drag for orbiting and scroll for zooming. Use keys to change quantum numbers or particle counts in real time.