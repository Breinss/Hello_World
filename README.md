# 3D Particle Simulation

This is my first ever C++ project—a 3D particle physics simulation!

---

## Project Overview

This application simulates particle physics in 3D space with the following features:

- **Real-time 3D rendering** with OpenGL
- **Particle-particle collision detection** and resolution
- **Interactive camera controls**
- **Adjustable physics parameters** (gravity, particle size, etc.)
- **Two boundary types**: box and sphere
- **Performance scaling** based on system capabilities

---

## First C++ Learning Project

This project represents my first foray into C++ programming. As a beginner, I learned:

- Object-oriented programming in C++
- Working with pointers and memory management
- OpenGL integration and shader programming
- Cross-platform development with CMake
- Physics simulation techniques

---

## Technical Features

### Physics Engine
- **Semi-implicit Euler integration**
- **Collision detection** with spatial partitioning
- **Fluid dynamics** using the SPH algorithm
- **Variable time stepping** for stability

### Graphics
- **Modern OpenGL rendering pipeline**
- **Custom GLSL shaders**
- **Point sprites** for efficient particle rendering
- **Real-time performance adaptation**

### User Interface
- **Interactive camera controls**
- **On-screen parameter display**
- **Keyboard shortcuts** for adjusting simulation parameters

---

## Dependencies

The project uses the following libraries and tools:

- **SFML 2.5**
- **OpenGL**
- **GLEW**
- **GLM**
- **ImGui & ImGui-SFML** (fetched by CMake)

---

## Building the Project

To build the project, ensure all dependencies are installed and follow the instructions in the `CMakeLists.txt` file.

---

## Controls

| Key         | Action                                   |
|-------------|-----------------------------------------|
| Arrow Keys  | Rotate camera                           |
| `R`         | Reset camera position                   |
| `P`         | Reset particles                         |
| `1` / `2`   | Decrease / Increase gravity             |
| `3` / `4`   | Decrease / Increase particle radius     |
| `5` / `6`   | Decrease / Increase restitution         |
| `7` / `8`   | Decrease / Increase boundary size       |
| `9` / `0`   | Decrease / Increase particle count      |
| `-` / `+`   | Decrease / Increase particle mass       |
| `F`         | Toggle friction                         |
| `L`         | Toggle fluid simulation                 |
| `B`         | Toggle boundary type (sphere/box)       |
| `T`         | Toggle performance mode                 |

---

## What I've Learned

As my first C++ project, this simulation taught me a tremendous amount about:

- C++ syntax and language features
- Building complex applications from scratch
- Graphics programming fundamentals
- Physics simulation algorithms
- Performance optimization techniques
- Modern C++ practices and tools

The particle simulation is both visually appealing and physically accurate, demonstrating principles of fluid dynamics and rigid body physics in an interactive environment.