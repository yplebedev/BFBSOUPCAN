// Port from wikipeda
// https://en.wikipedia.org/wiki/Halton_sequence#Implementation_in_pseudocode

#define PI 3.1515
#define PHI 1.6180339887

float halton(int b, int i) {
	float result = 0.;
	float f = 1.;
	while (i > 0)
	{
		f = f / float(b);
		result += f * float(i % b);
		i = i / b; 
        //index = int(floor(float(index) / float(base)));
	}
	return result;
}

uniform int random_value < source = "random"; min = 0; max = 4000; >;

// A BUNCH of code from UberRT, thanks yet again <3
int getPixelIndex(float2 coords, int2 size) {
    int2 pxCoords = coords*size;
    return pxCoords.x+pxCoords.y*size.x+random_value;
}

float randomValue(inout uint seed) {
    seed = seed * 747796405 + 2891336453;
    uint result = ((seed>>((seed>>28)+4))^seed)*277803737;
    result = (result>>22)^result;
    return result/4294967295.0;
}

float3 randomTriple(in out uint seed) {
	float3 v = 0;
	v.x = randomValue(seed);
	v.y = randomValue(seed);
	v.z = randomValue(seed);
	return v;
}

float boxmuller(float r0, float r1) {
	return 0.25 * sqrt(-log(r0 + 0.00001))*cos(PI * 2.0 * r1) + 0.5;
}

float3 randomGauss3(float2 coords, in out uint seed) {
	float3 r = randomTriple(seed);
	float ox = boxmuller(r.x, r.y);
	float oy = boxmuller(r.x, r.z);
	float oz = boxmuller(r.z, r.y);
	return 0.5 - 2 * float3(ox, oy, oz);
}

float3 randomVec3(float2 coords, in out uint seed) {
	while (true) {
		float3 rng = randomGauss3(coords, seed);
		if (length(rng) <= 1.0) return normalize(rng);
	}
}