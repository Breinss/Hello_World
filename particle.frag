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
