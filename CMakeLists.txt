cmake_minimum_required(VERSION 3.14)
project(ParticleSimulation VERSION 1.0.0 LANGUAGES CXX)

# C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find SFML packages
find_package(SFML 2.5 COMPONENTS graphics window system REQUIRED)

# Find OpenGL, GLUT, GLEW, and glm
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(GLEW REQUIRED)
find_package(glm REQUIRED) # Add glm

# FetchContent for ImGui and ImGui-SFML
include(FetchContent)

# Fetch ImGui first (using a compatible version)
FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG v1.86 # Using v1.86 which is compatible with ImGui-SFML v2.5
)

# Make ImGui available first
FetchContent_MakeAvailable(imgui)

# Now that ImGui is available, set the variable for ImGui-SFML
set(IMGUI_DIR ${imgui_SOURCE_DIR})

# Fetch ImGui-SFML
FetchContent_Declare(
    imgui_sfml
    GIT_REPOSITORY https://github.com/SFML/imgui-sfml.git
    GIT_TAG v2.5
)

# Configure ImGui-SFML
set(IMGUI_SFML_FIND_SFML OFF)
set(IMGUI_SFML_IMGUI_DEMO OFF) # Set to ON if you want the demo

# Make ImGui-SFML available
FetchContent_MakeAvailable(imgui_sfml)

# Define target executable
add_executable(simulation 
    main.cpp
    particle.cpp
)

# Include directories
target_include_directories(simulation PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${imgui_SOURCE_DIR}
)

# Link libraries
target_link_libraries(simulation PRIVATE
    ImGui-SFML::ImGui-SFML
    sfml-graphics
    sfml-window
    sfml-system
    OpenGL::GL
    GLUT::GLUT
    GLEW::GLEW
    glm::glm # Link glm
)