#version 410

in vec3 aPosition;
in vec3 aNormal;
in vec3 aColor;
in vec2 aTexCoord;
in vec3 aCentroid;
in vec3 aInstancedPosition;

uniform float uScale;
uniform mat4 uModelViewMatrix;
uniform mat4 uNormalMatrix;
uniform mat4 uProjectionMatrix;
uniform bool uUseTexture;
uniform sampler2D uTexture;
uniform bool uLight;

out vec3 vPosition;
out vec3 vNormal;
out vec3 vColor;
out vec2 vTexCoord;

void main() {
    vec4 modelViewPosition = uModelViewMatrix * vec4(uScale * (aPosition + aInstancedPosition) + (1 - uScale) * aCentroid, 1);
    vPosition = modelViewPosition.xyz;
    vNormal = normalize((uNormalMatrix * vec4(aNormal, 0)).xyz);
    vColor = aColor;
    vTexCoord = aTexCoord;
    gl_Position = uProjectionMatrix * modelViewPosition;
}
