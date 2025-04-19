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
