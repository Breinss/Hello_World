#ifndef PARTICLE_H
#define PARTICLE_H

class Particle {
private:
    double x, y, z;        // 3D Position
    double vx, vy, vz;     // 3D Velocity
    double ax, ay, az;     // 3D Acceleration
    double mass;           // Particle mass
    double radius;         // Particle radius
    double restitution;    // Coefficient of restitution (bounciness)
    double density;        // Fluid density at particle

public:
    Particle(double x_pos, double y_pos, double z_pos);
    
    // Getters
    double getX() const { return x; }
    double getY() const { return y; }
    double getZ() const { return z; }
    double getVX() const { return vx; }
    double getVY() const { return vy; }
    double getVZ() const { return vz; }
    double getRadius() const { return radius; }
    double getRestitution() const { return restitution; }
    double getMass() const { return mass; }
    double getDensity() const { return density; }
    bool isOnGround(double boundaryHeight, bool isSphereBoundary = false, double sphereRadius = 0.0) const;
    
    // Setters
    void setX(double new_x) { x = new_x; }
    void setY(double new_y) { y = new_y; }
    void setVX(double new_vx) { vx = new_vx; }
    void setVY(double new_vy) { vy = new_vy; }
    void setRadius(double new_radius) { radius = new_radius; }
    void setRestitution(double new_restitution) { restitution = new_restitution; }
    void setMass(double new_mass) { mass = new_mass; }
    void setDensity(double new_density) { density = new_density; }
    
    void applyForce(double fx, double fy, double fz);
    void applyGravity(double g);
    void update(double dt);
    void resetForces();
    void handleBoundaryCollision(double width, double height, double depth);
    void handleSphereBoundaryCollision(double radius);
    bool checkCollision(const Particle& other) const;
    void resolveCollision(Particle& other, bool fluidMode = false);
};

#endif