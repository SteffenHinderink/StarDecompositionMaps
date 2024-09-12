#version 410

in vec3 vPosition;
in vec3 vNormal;
in vec3 vColor;
in vec2 vTexCoord;

uniform float uScale;
uniform mat4 uModelViewMatrix;
uniform mat4 uNormalMatrix;
uniform mat4 uProjectionMatrix;
uniform bool uUseTexture;
uniform sampler2D uTexture;
uniform bool uLight;

out vec4 outColor;

void main() {
    vec3 color = vColor;
    if (uUseTexture) {
        color = texture(uTexture, vTexCoord).rgb;
    }
    if (uLight) {
        vec3 ca = 0.3 * color;
        vec3 cd = 0.7 * color;
        vec3 cs = 0.1 * vec3(1, 1, 1);
        float s = 10.0;
        vec3 n = normalize(vNormal);
        vec3 v = normalize(-vPosition);
        if (dot(n, v) < 0.0) {
            n = -n;
        }
        vec3 l = v;
        vec3 r = 2.0 * n * dot(n, l) - l;
        color = (ca + cd * max(0.0, dot(l, n)) + cs * pow(max(0.0, dot(v, r)), s)) * vec3(1, 1, 1);
    }
    outColor = vec4(color, 1);
}
