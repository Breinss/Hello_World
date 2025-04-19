#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#define GLEW_STATIC  // Use static linking for GLEW
#include <GL/glew.h> // GLEW must be included before other OpenGL headers
#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>
#include <GL/gl.h>
#include <GL/glu.h>
#include "particle.h"
#include <sstream>
#include <cmath>
#include <fstream> // For shader file reading
#include <glm/glm.hpp> // GLM for matrix operations
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <unistd.h> // For getcwd
#include <unordered_map>
#include <vector>
#include <tuple>
#include <iomanip>
// Shader loading utility function - prototype declaration
GLuint LoadShaders(const char * vertex_file_path, const char * fragment_file_path);
GLuint CreateDefaultShaders();

// Define a simple grid cell structure
struct GridCell {
    std::vector<size_t> particleIndices;
};

// Grid-based spatial partitioning - add before main()
class SpatialGrid {
private:
    std::unordered_map<std::string, GridCell> grid;
    float cellSize;

public:
    SpatialGrid(float size) : cellSize(size) {}

    std::string getCellKey(float x, float y, float z) {
        int i = static_cast<int>(floor(x / cellSize));
        int j = static_cast<int>(floor(y / cellSize));
        int k = static_cast<int>(floor(z / cellSize));
        return std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k);
    }

    void clear() {
        grid.clear();
    }

    void insertParticle(size_t index, float x, float y, float z) {
        std::string key = getCellKey(x, y, z);
        grid[key].particleIndices.push_back(index);
    }

    std::vector<size_t> getNeighborCandidates(float x, float y, float z) {
        std::vector<size_t> candidates;
        
        // Get particles in current cell and adjacent cells
        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
                for (int k = -1; k <= 1; k++) {
                    float nx = x + i * cellSize;
                    float ny = y + j * cellSize;
                    float nz = z + k * cellSize;
                    
                    std::string key = getCellKey(nx, ny, nz);
                    if (grid.find(key) != grid.end()) {
                        auto& cell = grid[key];
                        candidates.insert(candidates.end(), cell.particleIndices.begin(), cell.particleIndices.end());
                    }
                }
            }
        }
        
        return candidates;
    }
};

// Helper function to draw a wireframe cube
void drawBoundingBox(float size) {
    glLineWidth(1.0f);
    glColor3f(0.5f, 0.5f, 0.5f);  // Gray color for the wireframe
    
    glBegin(GL_LINES);
    
    // Bottom face
    glVertex3f(-size, -size, -size); glVertex3f(size, -size, -size);
    glVertex3f(size, -size, -size); glVertex3f(size, -size, size);
    glVertex3f(size, -size, size); glVertex3f(-size, -size, size);
    glVertex3f(-size, -size, size); glVertex3f(-size, -size, -size);
    
    // Top face
    glVertex3f(-size, size, -size); glVertex3f(size, size, -size);
    glVertex3f(size, size, -size); glVertex3f(size, size, size);
    glVertex3f(size, size, size); glVertex3f(-size, size, size);
    glVertex3f(-size, size, size); glVertex3f(-size, size, -size);
    
    // Connecting edges
    glVertex3f(-size, -size, -size); glVertex3f(-size, size, -size);
    glVertex3f(size, -size, -size); glVertex3f(size, size, -size);
    glVertex3f(size, -size, size); glVertex3f(size, size, size);
    glVertex3f(-size, -size, size); glVertex3f(-size, size, size);
    
    glEnd();
}

// Helper function to draw a wireframe sphere
void drawBoundingSphere(float radius, int segments = 32) {
    glLineWidth(1.0f);
    glColor3f(0.5f, 0.5f, 0.5f);  // Gray color for the wireframe
    
    // Draw longitude lines (vertical)
    for (int i = 0; i < segments; i++) {
        float angle = (2.0f * M_PI * i) / segments;
        float nextAngle = (2.0f * M_PI * (i + 1)) / segments;
        
        float x = cos(angle);
        float nextX = cos(nextAngle);
        float z = sin(angle);
        float nextZ = sin(nextAngle);
        
        glBegin(GL_LINE_LOOP);
        for (int j = 0; j < segments; j++) {
            float phi = (M_PI * j) / (segments - 1) - M_PI/2;
            float y = sin(phi);
            float r = cos(phi);
            
            glVertex3f(radius * r * x, radius * y, radius * r * z);
        }
        glEnd();
    }
    
    // Draw latitude lines (horizontal)
    for (int j = 1; j < segments-1; j++) {
        float phi = (M_PI * j) / (segments - 1) - M_PI/2;
        float y = sin(phi);
        float r = cos(phi);
        
        glBegin(GL_LINE_LOOP);
        for (int i = 0; i <= segments; i++) {
            float angle = (2.0f * M_PI * i) / segments;
            float x = cos(angle);
            float z = sin(angle);
            
            glVertex3f(radius * r * x, radius * y, radius * r * z);
        }
        glEnd();
    }
}

// Function to reset particles to random positions
void resetParticles(std::vector<Particle>& particles, float boundarySize, int count, bool isSphereBoundary) {
    particles.clear();
    for(int i = 0; i < count; i++) {
        float x, y, z;
        
        if (isSphereBoundary) {
            // For spherical boundary, distribute particles in the top hemisphere
            float theta = (rand() % 628) / 100.0f;  // 0 to 2π
            float phi = (rand() % 157) / 100.0f;     // 0 to π/2 (top hemisphere)
            float r = boundarySize * 0.8f * (rand() % 100) / 100.0f;  // Random radius within 80% of boundary
            
            x = r * cos(phi) * cos(theta);
            y = r * sin(phi) + boundarySize * 0.1f; // Keep particles in top portion
            z = r * cos(phi) * sin(theta);
        } else {
            // For box boundary, use the existing distribution at the top
            x = (rand() % 200 - 100) / 100.0f * boundarySize;
            y = boundarySize - (rand() % 20) / 100.0f * boundarySize;
            z = (rand() % 200 - 100) / 100.0f * boundarySize;
        }
        
        particles.push_back(Particle(x, y, z));
    }
}

int main(int argc, char** argv) {
    // Add these variables at the beginning of main()
    bool lowPerformanceMode = false;
    float performanceScale = 1.0f;
    int frameSkip = 0;
    int currentFrame = 0;
    sf::Clock performanceClock;

    // Create window with OpenGL context
    sf::ContextSettings settings;
    settings.depthBits = 24;
    settings.stencilBits = 8;
    settings.antialiasingLevel = 4;
    settings.majorVersion = 3; // Request OpenGL 3.3 context
    settings.minorVersion = 3;
    settings.attributeFlags = sf::ContextSettings::Default; // Use compatibility profile instead of Core

    sf::RenderWindow window(sf::VideoMode(800, 600), "3D Particle Simulation", sf::Style::Default, settings);
    window.setVerticalSyncEnabled(false);
    window.setFramerateLimit(60); // Cap at 60 FPS

    // Initialize GLEW after creating the window and OpenGL context
    glewExperimental = GL_TRUE; // Needed for core profile - set this BEFORE glewInit
    GLenum glewError = glewInit();
    if (glewError != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW: " << glewGetErrorString(glewError) << std::endl;
        return -1;
    }
    
    // Sometimes glewInit generates an OpenGL error - clear it
    glGetError();
    
    std::cout << "OpenGL Version: " << glGetString(GL_VERSION) << std::endl;
    std::cout << "GLSL Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
    std::cout << "Vendor: " << glGetString(GL_VENDOR) << std::endl;
    std::cout << "Renderer: " << glGetString(GL_RENDERER) << std::endl;

    // Load font for on-screen controls
    sf::Font font;
    if (!font.loadFromFile("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")) {
        // Try a fallback font if the first one fails
        if (!font.loadFromFile("/usr/share/fonts/truetype/freefont/FreeSans.ttf")) {
            std::cout << "Error loading font. Controls will not be displayed on screen." << std::endl;
        }
    }

    // Create text objects for on-screen controls
    sf::Text controlsText;
    controlsText.setFont(font);
    controlsText.setCharacterSize(16);
    controlsText.setFillColor(sf::Color::White);
    controlsText.setPosition(10, 10);
    
    sf::Text statsText;
    statsText.setFont(font);
    statsText.setCharacterSize(16);
    statsText.setFillColor(sf::Color::Yellow);
    statsText.setPosition(10, 300);

    // Set fluid parameters
    float smoothingRadius = 0.5f;  // Radius for fluid interaction
    float fluidStiffness = 100.0f;  // Pressure force coefficient
    float fluidRestDensity = 20.0f;  // Target density for fluid
    float fluidViscosity = 2.0f;  // How "thick" the fluid is
    bool fluidMode = true;  // Enable fluid simulation by default

    // Set up 3D perspective using GLM
    glEnable(GL_DEPTH_TEST);
    glm::mat4 ProjectionMatrix = glm::perspective(glm::radians(45.0f), (float)window.getSize().x / (float)window.getSize().y, 0.1f, 500.0f);

    // Seed random number generator
    srand(time(0));

    // Simulation parameters
    double dt = 0.008; // Balance between stability and speed
    double gravity = -9.81;
    float boundarySize = 5.0f;
    bool isSphereBoundary = true; // Default to sphere boundary
    
    // Particle properties - fluid particles work better with slightly larger sizes
    float particleRadius = 0.08f;
    float particleRestitution = 0.1f; // Lower restitution for fluid-like behavior
    float particleMass = 1.0f;
    int particleCount = 1000;

    // Add physics enhancement variables
    bool hasFriction = true;
    float frictionCoeff = 0.2f;

    // Create initial particles with appropriate boundary
    std::vector<Particle> particles;
    resetParticles(particles, boundarySize, particleCount, isSphereBoundary);

    // Initialize spatial grid with cell size slightly larger than particle diameter
    SpatialGrid spatialGrid(smoothingRadius * 2);

    // Static camera setup with rotation controls
    float cameraAngleX = 0.0f;
    float cameraAngleY = 0.0f;
    float rotationSpeed = 0.0f;
    float elevationSpeed = 0.0f;
    float rotationDamping = 0.95f; // Damping factor for smooth rotation

    // Keep the background color for better visibility
    glClearColor(0.0f, 0.0f, 0.1f, 1.0f);  // Dark blue background

    // Display help message
    std::cout << "\n=== PARTICLE SIMULATION CONTROLS ===\n";
    std::cout << "Arrow Keys: Rotate camera\n";
    std::cout << "R: Reset camera position\n";
    std::cout << "P: Reset particles\n";
    std::cout << "1/2: Decrease/Increase gravity\n";
    std::cout << "3/4: Decrease/Increase particle radius\n";
    std::cout << "5/6: Decrease/Increase restitution\n";
    std::cout << "7/8: Decrease/Increase boundary size\n";
    std::cout << "9/0: Decrease/Increase particle count\n";
    std::cout << "-/+: Decrease/Increase particle mass\n";
    std::cout << "F: Toggle friction (currently " << (hasFriction ? "ON" : "OFF") << ")\n";
    std::cout << "L: Toggle fluid simulation (currently " << (fluidMode ? "ON" : "OFF") << ")\n";
    std::cout << "B: Toggle boundary shape (currently " << (isSphereBoundary ? "SPHERE" : "BOX") << ")\n";
    std::cout << "T: Toggle performance mode\n";
    std::cout << "==================================\n\n";

    // --- GPU Instancing Setup ---
    std::cout << "Attempting to load shaders from current directory: " << std::endl;
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        std::cout << "Current working directory: " << cwd << std::endl;
    }
    
    std::cout << "Looking for: particle.vert and particle.frag" << std::endl;
    
    GLuint particleShaderProgram = LoadShaders("particle.vert", "particle.frag");
    if (particleShaderProgram == 0) {
        std::cerr << "Failed to load shaders. Using default fallback shaders." << std::endl;
        // Try fallback path - maybe they're in a shaders subdirectory?
        particleShaderProgram = LoadShaders("shaders/particle.vert", "shaders/particle.frag");
        
        if (particleShaderProgram == 0) {
            std::cerr << "Fallback shaders also failed. Creating default in-memory shaders." << std::endl;
            particleShaderProgram = CreateDefaultShaders();
            
            if (particleShaderProgram == 0) {
                std::cerr << "FATAL: Failed to create default shaders. Exiting." << std::endl;
                return -1;
            }
        }
    }

    // Get shader uniform locations
    GLuint projectionMatrixID = glGetUniformLocation(particleShaderProgram, "ProjectionMatrix");
    GLuint viewMatrixID = glGetUniformLocation(particleShaderProgram, "ViewMatrix");

    // --- VBOs for Instanced Data ---
    GLuint particlePositionBuffer;
    glGenBuffers(1, &particlePositionBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, particlePositionBuffer);
    glBufferData(GL_ARRAY_BUFFER, particleCount * sizeof(glm::vec4), NULL, GL_STREAM_DRAW); // Position (xyz) + Radius (w)

    GLuint particleColorBuffer;
    glGenBuffers(1, &particleColorBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, particleColorBuffer);
    glBufferData(GL_ARRAY_BUFFER, particleCount * sizeof(glm::vec4), NULL, GL_STREAM_DRAW); // Color (rgba)

    // --- VAO Setup ---
    GLuint particleVAO;
    glGenVertexArrays(1, &particleVAO);
    glBindVertexArray(particleVAO);

    // Attribute 0: Position + Radius (from particlePositionBuffer)
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, particlePositionBuffer);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glVertexAttribDivisor(0, 1); // Tell OpenGL this is an instanced vertex attribute.

    // Attribute 1: Color (from particleColorBuffer)
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, particleColorBuffer);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glVertexAttribDivisor(1, 1); // Tell OpenGL this is an instanced vertex attribute.

    glBindVertexArray(0); // Unbind VAO

    // Buffers for particle data
    std::vector<glm::vec4> particlePositions(particleCount);
    std::vector<glm::vec4> particleColors(particleCount);
    // --- End GPU Instancing Setup ---

    // Clock for timing
    sf::Clock deltaClock;
    sf::Clock fpsClock;
    float fpsUpdateInterval = 0.5f; // Update FPS in title every 0.5 seconds
    float currentFps = 0.0f;

    while (window.isOpen()) {
        sf::Time deltaTime = deltaClock.restart();
        float dt_seconds = deltaTime.asSeconds();
        
        // Calculate current FPS
        currentFps = 1.0f / dt_seconds;
        
        // Update window title with FPS periodically
        if (fpsClock.getElapsedTime().asSeconds() >= fpsUpdateInterval) {
            std::stringstream title;
            title << "3D Particle Simulation - FPS: " << std::fixed << std::setprecision(1) << currentFps;
            window.setTitle(title.str());
            fpsClock.restart();
        }

        // Add performance monitoring and quality adjustment
        if (performanceClock.getElapsedTime().asSeconds() >= 2.0) {
            float fps = 1.0f / dt_seconds;
            
            if (fps < 20.0f && !lowPerformanceMode) {
                // Switch to low performance mode if FPS drops
                lowPerformanceMode = true;
                particleRadius *= 1.5f;     // Larger particles for better visibility
                smoothingRadius *= 1.5f;    // Larger smoothing radius
                performanceScale = 0.5f;    // Half the particles
                frameSkip = 1;              // Skip physics calculations every other frame
                std::cout << "Switching to low performance mode" << std::endl;
            } else if (fps > 40.0f && lowPerformanceMode) {
                // Return to normal mode if performance improves
                lowPerformanceMode = false;
                particleRadius /= 1.5f;
                smoothingRadius /= 1.5f;
                performanceScale = 1.0f;
                frameSkip = 0;
                std::cout << "Returning to normal performance mode" << std::endl;
            }
            
            performanceClock.restart();
        }

        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
            
            // Handle key presses for camera rotation
            if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::Left)
                    rotationSpeed = 120.0f; // degrees per second
                if (event.key.code == sf::Keyboard::Right)
                    rotationSpeed = -120.0f;
                if (event.key.code == sf::Keyboard::Up)
                    elevationSpeed = 120.0f;
                if (event.key.code == sf::Keyboard::Down)
                    elevationSpeed = -120.0f;
                    
                // Reset camera rotation
                if (event.key.code == sf::Keyboard::R) {
                    cameraAngleX = 0.0f;
                    cameraAngleY = 0.0f;
                    rotationSpeed = 0.0f;
                    elevationSpeed = 0.0f;
                }
                
                // Reset particles with P key
                if (event.key.code == sf::Keyboard::P) {
                    resetParticles(particles, boundarySize, particleCount, isSphereBoundary);
                }
                
                // Parameter adjustment keys
                if (event.key.code == sf::Keyboard::Num1) gravity = std::max(-20.0, gravity - 0.5);
                if (event.key.code == sf::Keyboard::Num2) gravity = std::min(0.0, gravity + 0.5);
                
                if (event.key.code == sf::Keyboard::Num3) particleRadius = std::max(0.01f, particleRadius - 0.01f);
                if (event.key.code == sf::Keyboard::Num4) particleRadius = std::min(1.0f, particleRadius + 0.01f);
                
                if (event.key.code == sf::Keyboard::Num5) particleRestitution = std::max(0.0f, particleRestitution - 0.05f);
                if (event.key.code == sf::Keyboard::Num6) particleRestitution = std::min(1.0f, particleRestitution + 0.05f);
                
                if (event.key.code == sf::Keyboard::Num7) boundarySize = std::max(1.0f, boundarySize - 0.5f);
                if (event.key.code == sf::Keyboard::Num8) boundarySize = std::min(20.0f, boundarySize + 0.5f);
                
                if (event.key.code == sf::Keyboard::Num9) {
                    particleCount = std::max(100, particleCount - 100);
                    resetParticles(particles, boundarySize, particleCount, isSphereBoundary);
                }
                if (event.key.code == sf::Keyboard::Num0) {
                    particleCount = std::min(5000, particleCount + 100);
                    resetParticles(particles, boundarySize, particleCount, isSphereBoundary);
                }
                
                if (event.key.code == sf::Keyboard::Subtract) particleMass = std::max(0.1f, particleMass - 0.1f);
                if (event.key.code == sf::Keyboard::Add) particleMass = std::min(10.0f, particleMass + 0.1f);

                // Toggle friction
                if (event.key.code == sf::Keyboard::F) {
                    hasFriction = !hasFriction;
                    std::cout << "Friction " << (hasFriction ? "enabled" : "disabled") << std::endl;
                }

                // Toggle fluid mode
                if (event.key.code == sf::Keyboard::L) {
                    fluidMode = !fluidMode;
                    std::cout << "Fluid simulation " << (fluidMode ? "enabled" : "disabled") << std::endl;
                    if (fluidMode) {
                        particleRestitution = 0.1f; // Lower restitution for fluid
                    } else {
                        particleRestitution = 0.7f; // Higher for regular particles
                    }
                }

                // Toggle boundary type with B key
                if (event.key.code == sf::Keyboard::B) {
                    isSphereBoundary = !isSphereBoundary;
                    std::cout << "Boundary shape: " << (isSphereBoundary ? "SPHERE" : "BOX") << std::endl;
                }

                // Toggle performance mode with T key
                if (event.key.code == sf::Keyboard::T) {
                    lowPerformanceMode = !lowPerformanceMode;
                    if (lowPerformanceMode) {
                        particleCount /= 2;
                        resetParticles(particles, boundarySize, particleCount, isSphereBoundary);
                        std::cout << "Performance mode enabled: " << particleCount << " particles" << std::endl;
                    } else {
                        particleCount *= 2;
                        resetParticles(particles, boundarySize, particleCount, isSphereBoundary);
                        std::cout << "Performance mode disabled: " << particleCount << " particles" << std::endl;
                    }
                }
            }
            
            // Handle key releases for camera control
            if (event.type == sf::Event::KeyReleased) {
                if (event.key.code == sf::Keyboard::Left && rotationSpeed > 0)
                    rotationSpeed = 0.0f;
                if (event.key.code == sf::Keyboard::Right && rotationSpeed < 0)
                    rotationSpeed = 0.0f;
                if (event.key.code == sf::Keyboard::Up && elevationSpeed > 0)
                    elevationSpeed = 0.0f;
                if (event.key.code == sf::Keyboard::Down && elevationSpeed < 0)
                    elevationSpeed = 0.0f;
            }
        }

        // Apply damping and update camera angles
        rotationSpeed *= rotationDamping;
        elevationSpeed *= rotationDamping;
        
        cameraAngleY += rotationSpeed * dt_seconds;
        cameraAngleX += elevationSpeed * dt_seconds;
        
        // Limit vertical rotation to avoid flipping
        if (cameraAngleX > 89.0f) cameraAngleX = 89.0f;
        if (cameraAngleX < -89.0f) cameraAngleX = -89.0f;

        // Skip physics on some frames
        currentFrame++;
        if (frameSkip > 0 && currentFrame % (frameSkip + 1) != 0) {
            // Skip physics calculations but still render
        } else {
            // Update spatial grid at the beginning of physics calculations
            spatialGrid.clear();
            for (size_t i = 0; i < particles.size(); i++) {
                spatialGrid.insertParticle(i, particles[i].getX(), particles[i].getY(), particles[i].getZ());
            }

            // Update particle properties - now with more realistic physics
            for (auto& p : particles) {
                p.resetForces();
                p.applyGravity(gravity);
                
                // Apply friction if enabled
                if (hasFriction && p.isOnGround(boundarySize)) {
                    double vx = p.getVX();
                    double vz = p.getVZ();
                    double speed = std::sqrt(vx*vx + vz*vz);
                    
                    if (speed > 0.001) {
                        double frictionFx = -vx/speed * frictionCoeff * 9.81 * p.getMass();
                        double frictionFz = -vz/speed * frictionCoeff * 9.81 * p.getMass();
                        p.applyForce(frictionFx, 0, frictionFz);
                    }
                }
            }

            // Calculate fluid forces if fluid mode is enabled
            if (fluidMode) {
                // First pass: calculate density using spatial grid
                for(size_t i = 0; i < particles.size(); i++) {
                    double density = 0.0;
                    auto candidates = spatialGrid.getNeighborCandidates(
                        particles[i].getX(), particles[i].getY(), particles[i].getZ());
                    
                    for (auto j : candidates) {
                        if (i == j) continue;
                        
                        // Calculate distance between particles
                        double dx = particles[j].getX() - particles[i].getX();
                        double dy = particles[j].getY() - particles[i].getY();
                        double dz = particles[j].getZ() - particles[i].getZ();
                        double distSquared = dx*dx + dy*dy + dz*dz;
                        
                        // Skip square root calculation when possible
                        if (distSquared < smoothingRadius * smoothingRadius) {
                            double dist = std::sqrt(distSquared);
                            double influence = 1.0 - dist/smoothingRadius;
                            density += influence * influence * particles[j].getMass();
                        }
                    }
                    
                    particles[i].setDensity(density);
                }
                
                // Second pass: apply fluid forces using spatial grid
                for(size_t i = 0; i < particles.size(); i++) {
                    auto candidates = spatialGrid.getNeighborCandidates(
                        particles[i].getX(), particles[i].getY(), particles[i].getZ());
                    
                    for(auto j : candidates) {
                        if (i >= j) continue; // Process each pair only once
                        
                        // Calculate distance
                        double dx = particles[j].getX() - particles[i].getX();
                        double dy = particles[j].getY() - particles[i].getY();
                        double dz = particles[j].getZ() - particles[i].getZ();
                        double distSquared = dx*dx + dy*dy + dz*dz;
                        
                        if (distSquared >= smoothingRadius * smoothingRadius || distSquared < 0.0001) continue;
                        
                        double dist = std::sqrt(distSquared);
                        
                        // Rest of fluid calculation remains the same...
                        double nx = dx / dist;
                        double ny = dy / dist;
                        double nz = dz / dist;
                        
                        double pressureI = fluidStiffness * (particles[i].getDensity() - fluidRestDensity);
                        double pressureJ = fluidStiffness * (particles[j].getDensity() - fluidRestDensity);
                        double pressureMagnitude = (pressureI + pressureJ) * 0.5f;
                        
                        double dvx = particles[j].getVX() - particles[i].getVX();
                        double dvy = particles[j].getVY() - particles[i].getVY();
                        double dvz = particles[j].getVZ() - particles[i].getVZ();
                        
                        double influence = 1.0 - dist/smoothingRadius;
                        influence = influence * influence;
                        
                        double fx = (nx * pressureMagnitude - dvx * fluidViscosity) * influence;
                        double fy = (ny * pressureMagnitude - dvy * fluidViscosity) * influence;
                        double fz = (nz * pressureMagnitude - dvz * fluidViscosity) * influence;
                        
                        particles[i].applyForce(-fx, -fy, -fz);
                        particles[j].applyForce(fx, fy, fz);
                    }
                }
            }

            // Update all particles based on calculated forces
            for(size_t i = 0; i < particles.size(); i++) {
                particles[i].setRadius(particleRadius);
                particles[i].setRestitution(particleRestitution);
                particles[i].setMass(particleMass);
                particles[i].update(dt);
                
                // Handle boundary collisions with appropriate type
                if (isSphereBoundary) {
                    particles[i].handleSphereBoundaryCollision(boundarySize);
                } else {
                    particles[i].handleBoundaryCollision(boundarySize, boundarySize, boundarySize);
                }
            }

            // Check collisions using the grid
            for (size_t i = 0; i < particles.size(); i++) {
                auto candidates = spatialGrid.getNeighborCandidates(
                    particles[i].getX(), particles[i].getY(), particles[i].getZ());
                
                for (auto j : candidates) {
                    if (i != j && i < j) { // Avoid duplicate checks
                        if (particles[i].checkCollision(particles[j])) {
                            particles[i].resolveCollision(particles[j], fluidMode);
                        }
                    }
                }
            }
        }
        
        // --- Update Particle Data for GPU ---
        if (particles.size() != particlePositions.size()) {
             // Resize buffers if particle count changed
            particlePositions.resize(particles.size());
            particleColors.resize(particles.size());

            glBindBuffer(GL_ARRAY_BUFFER, particlePositionBuffer);
            glBufferData(GL_ARRAY_BUFFER, particles.size() * sizeof(glm::vec4), NULL, GL_STREAM_DRAW);
            glBindBuffer(GL_ARRAY_BUFFER, particleColorBuffer);
            glBufferData(GL_ARRAY_BUFFER, particles.size() * sizeof(glm::vec4), NULL, GL_STREAM_DRAW);
        }

        for(size_t i = 0; i < particles.size(); ++i) {
            particlePositions[i] = glm::vec4(particles[i].getX(), particles[i].getY(), particles[i].getZ(), particles[i].getRadius());

            if (fluidMode) {
                float densityFactor = std::min(1.0f, (float)(particles[i].getDensity() / (fluidRestDensity * 2.0f)));
                particleColors[i] = glm::vec4(0.0f, densityFactor, 1.0f, 1.0f); // Blue based on density
            } else {
                particleColors[i] = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f); // White
            }
        }

        // Update VBOs
        glBindBuffer(GL_ARRAY_BUFFER, particlePositionBuffer);
        glBufferSubData(GL_ARRAY_BUFFER, 0, particles.size() * sizeof(glm::vec4), particlePositions.data());

        glBindBuffer(GL_ARRAY_BUFFER, particleColorBuffer);
        glBufferSubData(GL_ARRAY_BUFFER, 0, particles.size() * sizeof(glm::vec4), particleColors.data());

        glBindBuffer(GL_ARRAY_BUFFER, 0); // Unbind buffer
        // --- End Update Particle Data ---

        // Update stats text
        std::stringstream stats;
        stats << "Particles: " << particles.size() << "\n"
              << "Gravity: " << gravity << "\n"
              << "Radius: " << particleRadius << "\n"
              << "Restitution: " << particleRestitution << "\n"
              << "Mass: " << particleMass << "\n"
              << "Friction: " << (hasFriction ? "ON" : "OFF") << "\n"
              << "Fluid Mode: " << (fluidMode ? "ON" : "OFF") << "\n"
              << "Boundary: " << (isSphereBoundary ? "SPHERE" : "BOX");
        statsText.setString(stats.str());

        // Build controls info text
        std::stringstream controls;
        controls << "=== CONTROLS ===\n"
                 << "Arrow Keys: Rotate camera\n"
                 << "R: Reset camera\n"
                 << "P: Reset particles\n"
                 << "1/2: Adjust gravity\n"
                 << "3/4: Adjust radius\n"
                 << "5/6: Adjust restitution\n"
                 << "7/8: Adjust boundary\n"
                 << "9/0: Adjust count\n"
                 << "-/+: Adjust mass\n"
                 << "F: Toggle friction\n"
                 << "L: Toggle fluid mode\n"
                 << "B: Toggle boundary\n"
                 << "T: Toggle performance mode";
        controlsText.setString(controls.str());

        // Start fresh frame
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Camera setup with GLM
        glm::mat4 ViewMatrix = glm::mat4(1.0f); // Identity
        
        // Make camera back up based on boundary size
        float cameraDistance = boundarySize * 4.0f; // looks about right
        ViewMatrix = glm::translate(ViewMatrix, glm::vec3(0.0f, 0.0f, -cameraDistance)); 
        
        ViewMatrix = glm::rotate(ViewMatrix, glm::radians(cameraAngleX), glm::vec3(1.0f, 0.0f, 0.0f)); // x rotation
        ViewMatrix = glm::rotate(ViewMatrix, glm::radians(cameraAngleY), glm::vec3(0.0f, 1.0f, 0.0f)); // y rotation
        // Camera done

        // Draw the boundary box/sphere & coordinate axes
        // Need to switch from shaders back to old-school GL

        // Properly transition from shader-based to fixed-function rendering
        glUseProgram(0);
        glBindVertexArray(0);
        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);
        
        // Push matrices to save state before changing them
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadMatrixf(glm::value_ptr(ProjectionMatrix));
        
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadMatrixf(glm::value_ptr(ViewMatrix));
        
        // Enable depth testing for proper occlusion
        glEnable(GL_DEPTH_TEST);
        
        // Draw appropriate boundary
        if (isSphereBoundary) {
            drawBoundingSphere(boundarySize);
        } else {
            drawBoundingBox(boundarySize);
        }

        // Keep the coordinate axes for better orientation
        glLineWidth(3.0f);  // Make lines thicker
        glBegin(GL_LINES);
        // X-axis (red)
        glColor3f(1.0f, 0.0f, 0.0f);
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(boundarySize, 0.0f, 0.0f);
        // Y-axis (green)
        glColor3f(0.0f, 1.0f, 0.0f);
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(0.0f, boundarySize, 0.0f);
        // Z-axis (blue)
        glColor3f(0.0f, 0.0f, 1.0f);
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(0.0f, 0.0f, boundarySize);
        glEnd();

        // Keep the origin marker
        glPointSize(5.0f);
        glBegin(GL_POINTS);
        glColor3f(1.0f, 0.0f, 0.0f);  // Red color for origin
        glVertex3f(0.0f, 0.0f, 0.0f);
        glEnd();
        
        // Restore matrices after fixed-function rendering
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        // --- End Boundary and Axes ---

        // --- Render Particles using Instancing ---
        glUseProgram(particleShaderProgram);

        // Send MVP matrices to the shader
        glUniformMatrix4fv(projectionMatrixID, 1, GL_FALSE, glm::value_ptr(ProjectionMatrix));
        glUniformMatrix4fv(viewMatrixID, 1, GL_FALSE, glm::value_ptr(ViewMatrix));

        glBindVertexArray(particleVAO);
        glDrawArraysInstanced(GL_POINTS, 0, 1, particles.size()); // Draw 1 point per instance

        glBindVertexArray(0); // Unbind VAO
        glUseProgram(0); // Unbind shader program
        // --- End Particle Rendering ---

        // --- Render UI Text ---
        // Clear any errors that might have occurred before calling SFML functions
        glGetError();
        
        // Switch to 2D mode for text rendering with more careful state management
        window.pushGLStates();
        
        // Explicitly disable features that might interfere with SFML's rendering
        glDisable(GL_DEPTH_TEST);
        glActiveTexture(GL_TEXTURE0); // Reset active texture unit
        
        // Draw controls and stats
        window.draw(controlsText);
        window.draw(statsText);

        // Return to 3D mode
        window.popGLStates();
        
        // Restore OpenGL state needed for our 3D rendering
        glEnable(GL_DEPTH_TEST);
        // --- End UI Text ---

        // Display frame
        window.display();
    }

    // Cleanup GPU resources
    glDeleteBuffers(1, &particlePositionBuffer);
    glDeleteBuffers(1, &particleColorBuffer);
    glDeleteVertexArrays(1, &particleVAO);
    glDeleteProgram(particleShaderProgram);

    return 0;
}

// Default shader creation function
GLuint CreateDefaultShaders() {
    std::cout << "Creating default shaders in memory..." << std::endl;

    // Simple vertex shader for particle rendering
    const char* vertexShaderSource = R"glsl(
        #version 330 core
        layout(location = 0) in vec4 position_size; // xyz = position, w = size
        layout(location = 1) in vec4 color;
        
        uniform mat4 ProjectionMatrix;
        uniform mat4 ViewMatrix;
        
        out vec4 particleColor;
        
        void main() {
            particleColor = color;
            gl_Position = ProjectionMatrix * ViewMatrix * vec4(position_size.xyz, 1.0);
            gl_PointSize = position_size.w * 100.0; // Scale up for visibility
        }
    )glsl";

    // Simple fragment shader for particle rendering
    const char* fragmentShaderSource = R"glsl(
        #version 330 core
        in vec4 particleColor;
        out vec4 FragColor;
        
        void main() {
            // Create a circular particle
            vec2 circCoord = 2.0 * gl_PointCoord - 1.0;
            if (dot(circCoord, circCoord) > 1.0) {
                discard;
            }
            FragColor = particleColor;
        }
    )glsl";

    // Create and compile vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);
    
    // Check vertex shader compilation
    GLint success;
    GLchar infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cerr << "ERROR: Vertex shader compilation failed\n" << infoLog << std::endl;
        return 0;
    }

    // Create and compile fragment shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);
    
    // Check fragment shader compilation
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        std::cerr << "ERROR: Fragment shader compilation failed\n" << infoLog << std::endl;
        return 0;
    }

    // Create shader program and link shaders
    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    
    // Check linking
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cerr << "ERROR: Shader program linking failed\n" << infoLog << std::endl;
        return 0;
    }

    // Clean up
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    
    std::cout << "Default shaders created successfully." << std::endl;
    return shaderProgram;
}

// Enhanced Shader loading utility function - implementation after main
GLuint LoadShaders(const char * vertex_file_path, const char * fragment_file_path){
    std::cout << "Attempting to load shaders from: " << std::endl;
    std::cout << " - Vertex: " << vertex_file_path << std::endl;
    std::cout << " - Fragment: " << fragment_file_path << std::endl;

    // Create the shaders
    GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
    GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);

    // Read the Vertex Shader code from the file
    std::string VertexShaderCode;
    std::ifstream VertexShaderStream(vertex_file_path, std::ios::in);
    if(VertexShaderStream.is_open()){
        std::stringstream sstr;
        sstr << VertexShaderStream.rdbuf();
        VertexShaderCode = sstr.str();
        VertexShaderStream.close();
        std::cout << "Successfully opened vertex shader file" << std::endl;
    }else{
        std::cerr << "Error: Impossible to open " << vertex_file_path << ". Are you in the right directory?" << std::endl;
        return 0;
    }

    // Read the Fragment Shader code from the file
    std::string FragmentShaderCode;
    std::ifstream FragmentShaderStream(fragment_file_path, std::ios::in);
    if(FragmentShaderStream.is_open()){
        std::stringstream sstr;
        sstr << FragmentShaderStream.rdbuf();
        FragmentShaderCode = sstr.str();
        FragmentShaderStream.close();
        std::cout << "Successfully opened fragment shader file" << std::endl;
    } else {
        std::cerr << "Error: Impossible to open " << fragment_file_path << ". Are you in the right directory?" << std::endl;
        return 0;
    }

    GLint Result = GL_FALSE;
    int InfoLogLength;

    // Compile Vertex Shader
    std::cout << "Compiling shader: " << vertex_file_path << std::endl;
    char const * VertexSourcePointer = VertexShaderCode.c_str();
    glShaderSource(VertexShaderID, 1, &VertexSourcePointer , NULL);
    glCompileShader(VertexShaderID);

    // Check Vertex Shader
    glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    if (InfoLogLength > 0){
        std::vector<char> VertexShaderErrorMessage(InfoLogLength+1);
        glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
        std::cerr << "Vertex shader compilation error: " << &VertexShaderErrorMessage[0] << std::endl;
    }

    // Compile Fragment Shader
    std::cout << "Compiling shader: " << fragment_file_path << std::endl;
    char const * FragmentSourcePointer = FragmentShaderCode.c_str();
    glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer , NULL);
    glCompileShader(FragmentShaderID);

    // Check Fragment Shader
    glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    if (InfoLogLength > 0){
        std::vector<char> FragmentShaderErrorMessage(InfoLogLength+1);
        glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
        std::cerr << "Fragment shader compilation error: " << &FragmentShaderErrorMessage[0] << std::endl;
    }

    // Link the program
    std::cout << "Linking shader program" << std::endl;
    GLuint ProgramID = glCreateProgram();
    glAttachShader(ProgramID, VertexShaderID);
    glAttachShader(ProgramID, FragmentShaderID);
    glLinkProgram(ProgramID);

    // Check the program
    glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
    glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    if (InfoLogLength > 0){
        std::vector<char> ProgramErrorMessage(InfoLogLength+1);
        glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
        std::cerr << "Shader program linking error: " << &ProgramErrorMessage[0] << std::endl;
    }

    if (Result != GL_TRUE) {
        std::cerr << "ERROR: Shader program linking failed!" << std::endl;
        return 0;
    }

    std::cout << "Shader program created successfully with ID: " << ProgramID << std::endl;

    glDetachShader(ProgramID, VertexShaderID);
    glDetachShader(ProgramID, FragmentShaderID);

    glDeleteShader(VertexShaderID);
    glDeleteShader(FragmentShaderID);

    return ProgramID;
}