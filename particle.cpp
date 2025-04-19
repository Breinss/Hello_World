#include "particle.h"
#include <cstdlib>
#include <cmath>

Particle::Particle(double x_pos, double y_pos, double z_pos) {
    x = x_pos;
    y = y_pos;
    z = z_pos;
    vx = (rand() % 10 - 5) / 100.0;
    vy = 0.0; // No initial vertical velocity
    vz = (rand() % 10 - 5) / 100.0;
    ax = ay = az = 0.0;
    mass = 1.0;
    radius = 0.08;  // Slightly larger default for fluid behavior
    restitution = 0.1;  // Low restitution for fluid-like behavior
    density = 0.0;
}

void Particle::applyForce(double fx, double fy, double fz) {
    ax += fx / mass;
    ay += fy / mass;
    az += fz / mass;
}

void Particle::applyGravity(double g) {
    ay += g;
}

void Particle::update(double dt) {
    // Semi-implicit Euler integration (more stable)
    vx += ax * dt;
    vy += ay * dt;
    vz += az * dt;
    
    // Limit maximum velocity for stability
    double maxVelocity = 20.0;
    double speedSquared = vx*vx + vy*vy + vz*vz;
    if (speedSquared > maxVelocity*maxVelocity) {
        double speed = std::sqrt(speedSquared);
        double scale = maxVelocity / speed;
        vx *= scale;
        vy *= scale;
        vz *= scale;
    }
    
    // Update position using velocity
    x += vx * dt;
    y += vy * dt;
    z += vz * dt;
}

void Particle::resetForces() {
    ax = ay = az = 0.0;
}

void Particle::handleBoundaryCollision(double width, double height, double depth) {
    // Right wall
    if (x + radius > width) {
        x = width - radius;
        vx = -vx * restitution;
    }
    // Left wall
    if (x - radius < -width) {
        x = -width + radius;
        vx = -vx * restitution;
    }
    // Ceiling
    if (y + radius > height) {
        y = height - radius;
        vy = -vy * restitution;
    }
    // Floor - enhanced floor collision
    if (y - radius < -height) {
        y = -height + radius;
        
        // Only apply extra damping when hitting the floor at speed
        if (fabs(vy) > 0.5) {
            vy = -vy * restitution * 0.95; // Extra damping on fast impacts
        } else {
            vy = -vy * restitution;
            
            // Apply extra horizontal damping when nearly at rest
            if (fabs(vy) < 0.1) {
                vx *= 0.95; // More damping
                vz *= 0.95;
                
                // Set to true rest if very slow
                if (fabs(vy) < 0.01 && sqrt(vx*vx + vz*vz) < 0.01) {
                    vy = 0;
                }
            }
        }
    }
    // Front wall
    if (z + radius > depth) {
        z = depth - radius;
        vz = -vz * restitution;
    }
    // Back wall
    if (z - radius < -depth) {
        z = -depth + radius;
        vz = -vz * restitution;
    }
}

// New method to handle collisions with a spherical boundary
void Particle::handleSphereBoundaryCollision(double boundaryRadius) {
    // Calculate distance from origin (center of sphere)
    double distSquared = x*x + y*y + z*z;
    double dist = std::sqrt(distSquared);
    
    // Check if particle is outside or touching the boundary
    if (dist + radius > boundaryRadius) {
        // Calculate normal vector from center to particle
        double nx = x / dist;
        double ny = y / dist;
        double nz = z / dist;
        
        // Move particle back to boundary
        double penetration = dist + radius - boundaryRadius;
        x -= nx * penetration;
        y -= ny * penetration;
        z -= nz * penetration;
        
        // Calculate velocity along normal
        double velAlongNormal = vx * nx + vy * ny + vz * nz;
        
        // Reflect velocity
        vx = vx - 2 * velAlongNormal * nx;
        vy = vy - 2 * velAlongNormal * ny;
        vz = vz - 2 * velAlongNormal * nz;
        
        // Apply restitution
        vx *= restitution;
        vy *= restitution;
        vz *= restitution;
        
        // Check if particle is near the bottom of sphere for floor friction
        if (y < -boundaryRadius * 0.9) {
            vx *= 0.95; // Floor friction
            vz *= 0.95;
        }
    }
}

bool Particle::checkCollision(const Particle& other) const {
    double dx = other.x - x;
    double dy = other.y - y;
    double dz = other.z - z;
    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
    return distance < (radius + other.radius);
}

void Particle::resolveCollision(Particle& other, bool fluidMode) {
    double dx = other.x - x;
    double dy = other.y - y;
    double dz = other.z - z;
    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    if (distance < (radius + other.radius)) {
        // Normalize collision vector
        double nx = dx / distance;
        double ny = dy / distance;
        double nz = dz / distance;
        
        // Separation distance to resolve overlap
        double overlap = radius + other.radius - distance;
        
        // Minimal jitter for stability
        double jitter = (rand() % 2) / 1000.0;
        overlap *= (1.0 + jitter);
        
        // In fluid mode, allow some overlap for softer collisions
        double correctionFactor = fluidMode ? 0.5 : 1.0;
        
        // Move particles apart based on their masses
        double totalMass = mass + other.mass;
        double ratio1 = other.mass / totalMass;
        double ratio2 = mass / totalMass;
        
        double moveX = (overlap * nx) / 2.0 * correctionFactor;
        double moveY = (overlap * ny) / 2.0 * correctionFactor;
        double moveZ = (overlap * nz) / 2.0 * correctionFactor;
        
        x -= moveX * ratio1;
        y -= moveY * ratio1;
        z -= moveZ * ratio1;
        other.x += moveX * ratio2;
        other.y += moveY * ratio2;
        other.z += moveZ * ratio2;
        
        // Calculate relative velocity
        double relativeVx = other.vx - vx;
        double relativeVy = other.vy - vy;
        double relativeVz = other.vz - vz;
        
        // Calculate relative velocity along normal
        double velAlongNormal = relativeVx * nx + relativeVy * ny + relativeVz * nz;
        
        // Don't resolve if particles are moving apart
        if (velAlongNormal > 0) return;
        
        // For fluid-like behavior, increase damping but allow transfer of momentum
        double effectiveRestitution = fluidMode ? restitution * 0.5 : restitution;
        
        // Calculate impulse
        double j = -(1 + effectiveRestitution) * velAlongNormal;
        j /= 1/mass + 1/other.mass;
        
        // Apply impulse
        double impulseX = j * nx;
        double impulseY = j * ny;
        double impulseZ = j * nz;
        
        // In fluid mode, add tangential velocity transfer
        if (fluidMode) {
            // Calculate tangential component
            double dotProduct = relativeVx * nx + relativeVy * ny + relativeVz * nz;
            double tangentialVx = relativeVx - nx * dotProduct;
            double tangentialVy = relativeVy - ny * dotProduct;
            double tangentialVz = relativeVz - nz * dotProduct;
            
            // Add percentage of tangential velocity
            double tangentialFactor = 0.4; // Amount of sideways momentum transfer
            impulseX += tangentialVx * tangentialFactor;
            impulseY += tangentialVy * tangentialFactor;
            impulseZ += tangentialVz * tangentialFactor;
        }
        
        vx -= impulseX / mass;
        vy -= impulseY / mass;
        vz -= impulseZ / mass;
        other.vx += impulseX / other.mass;
        other.vy += impulseY / other.mass;
        other.vz += impulseZ / other.mass;
    }
}

bool Particle::isOnGround(double boundaryHeight, bool isSphereBoundary, double sphereRadius) const {
    if (isSphereBoundary) {
        // In sphere, "ground" is the lower hemisphere
        double distSquared = x*x + y*y + z*z;
        double dist = std::sqrt(distSquared);
        return (dist + radius > sphereRadius * 0.95) && (y < -sphereRadius * 0.7);
    } else {
        // In box, ground is the bottom face
        return y - radius <= -boundaryHeight + 0.01;
    }
}