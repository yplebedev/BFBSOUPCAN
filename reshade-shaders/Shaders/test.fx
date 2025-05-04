// VBGTAO
// MIT License...
/* Copyright (c)2025 Yaraslau Lebedzeu.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#include "ReShade.fxh"

#ifndef PI
	#define PI 3.14159265358979
#endif

#define getDepth ReShade::GetLinearizedDepth

uniform int framecount < source = "framecount"; >;

float3 getWorldPosition(float2 texcoord) {
	float depth = getDepth(texcoord);
    float3 result = float3((texcoord - 0.5) * depth, depth);
    result *= float3(BUFFER_WIDTH, BUFFER_HEIGHT, RESHADE_DEPTH_LINEARIZATION_FAR_PLANE);
    return result;
}

// It's slower because it samples around a lot for each pixel.
// Render "nicer" normals to a separate RT and fetch it later
float3 getWorldNormal(float2 texcoord) {
    float3 offset = float3(BUFFER_PIXEL_SIZE, 0);
    float2 posCenter = texcoord.xy;
    float2 posNorth  = posCenter - offset.zy;
    float2 posEast   = posCenter + offset.xz;
	float2 posSouth  = posCenter + offset.zy;
	float2 posWest   = posCenter - offset.xz;
	float  depthCenter = getDepth(posCenter);
	float  depthNorth = getDepth(posNorth);
	float  depthEast = getDepth(posEast);
	float  depthSouth = getDepth(posSouth);
	float  depthWest = getDepth(posWest);
	bool2 dir = bool2(abs(depthCenter - depthNorth) < abs(depthCenter - depthSouth), abs(depthCenter - depthEast) < abs(depthCenter - depthWest));
	
	float3 vertCenter = getWorldPosition(posCenter);

	if (dir.x && dir.y) {
		float3 vertNorth  = getWorldPosition(posNorth);
    	float3 vertEast   = getWorldPosition(posEast);
    	return -normalize(cross(vertCenter - vertNorth, vertCenter - vertEast));
	} else if (!dir.x && !dir.y) {
		float3 vertSouth  = getWorldPosition(posSouth);
    	float3 vertWest   = getWorldPosition(posWest);
    	return -normalize(cross(vertCenter - vertSouth, vertCenter - vertWest));		
	} else if (dir.x) {
		float3 vertNorth  = getWorldPosition(posNorth);
		float3 vertWest   = getWorldPosition(posWest);
		return -normalize(-cross(vertCenter - vertNorth, vertCenter - vertWest));
	} else {
		float3 vertSouth  = getWorldPosition(posSouth);
		float3 vertEast   = getWorldPosition(posEast);
		return -normalize(-cross(vertCenter - vertSouth, vertCenter - vertEast));
	}
}

// Number of set bits (slower than hardware implementation).
// Replace with a fast hardware call.
uint countbits(uint v) {
	uint c; 
	c = v - ((v >> 1) & 0x55555555);
	c = ((c >> 2) & 0x33333333) + (c & 0x33333333);
	c = ((c >> 4) + c) & 0x0F0F0F0F;
	c = ((c >> 8) + c) & 0x00FF00FF;
	c = ((c >> 16) + c) & 0x0000FFFF;
	return c;
}

#define SLICES 2

#define STEPS 5

#define SECTORS 32

#define R 800.0 // dear lawd almighty 

#define R_MAX_CLAMP 1200

#define FAR_CLIP (RESHADE_DEPTH_LINEARIZATION_FAR_PLANE-1)

#define THICKNESS 0.3

// Modified from https://www.shadertoy.com/view/4djSRW
// It's called once per pixel.
float2 fastHash(float2 p) {
	float3 p3 = float3(p * BUFFER_SCREEN_SIZE, framecount & 0xFFFF);
	p3 = frac(p3 * float3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx + 33.33);
    return frac((p3.xy + p3.yz) * p3.zx);
}

uint sliceSteps(float3 positionVS, float3 V, float2 start, float2 rayDir, float t, float step, float samplingDirection, float N, uint bitfield) {
    for (uint i = 0; i < STEPS; i++, t += step) {
        float2 samplePos = clamp(start + t * rayDir, 1, BUFFER_SCREEN_SIZE - 2);
        float3 samplePosVS = getWorldPosition(samplePos.xy / BUFFER_SCREEN_SIZE);
        float3 delta = samplePosVS - positionVS;
	
	    float2 fb = acos(float2(dot(normalize(delta), V), dot(normalize(delta - V * THICKNESS), V)));
	    fb = saturate(((samplingDirection * -fb) - N + PI/2) / PI);
	    fb = samplingDirection >= 0 ? fb.yx : fb.xy;

   		uint a = fb.x * SECTORS;
    	uint b = ceil((fb.y - fb.x) * SECTORS);
    	bitfield |= b > 0 ? (0xFFFFFFFF >> (SECTORS - b)) << a : 0;
    }
    return bitfield;
}

float gtao(float2 uv, float2 vpos) {
	float2 random = fastHash(uv);

	float2 start = uv * BUFFER_SCREEN_SIZE; // this is vpos
	float3 positionVS = getWorldPosition(uv);
	positionVS.z *= 0.9999; // Move center pixel towards camera a bit.

	float ao = 0.0;
	
	float3 V = normalize(-positionVS);
	float3 normalVS = getWorldNormal(uv);

    float step = max(1.0, clamp(R / positionVS.z, STEPS, R_MAX_CLAMP) / (STEPS + 1));

	for(float slice = 0.0; slice < 1.0; slice += 1.0 / SLICES) {
		float phi = PI * frac(slice + random.x);
		float2 direction = float2(cos(phi), sin(phi));

		float3 sliceN = normalize(cross(float3(direction, 0.0), V));
		float3 projN = normalVS - sliceN * dot(normalVS, sliceN);
		float cosN = dot(normalize(projN), V);

		float N = -sign(dot(projN, cross(V, sliceN))) * acos(cosN);
		
		uint aoBF = 0;
		float offset = max(random.y * step, length(BUFFER_PIXEL_SIZE));
		aoBF = sliceSteps(positionVS, V, start, direction, offset, step, 1, N, aoBF);
		aoBF = sliceSteps(positionVS, V, start, -direction, offset, step, -1, N, aoBF);

		ao += 1.0 - float(countbits(aoBF)) / float(SECTORS);
	}
	ao /= SLICES;
	return positionVS.z > FAR_CLIP || ao < -0.001 ? 1.0 : ao;
}

float4 main(float4 vpos : SV_Position, float2 uv : TEXCOORD) : SV_Target {
	float ao = gtao(uv, vpos.xy);
	return float4(ao.xxx, 1);
}

technique Test {
	pass Main { 
		VertexShader = PostProcessVS;
		PixelShader = main;
	}
}
