/// My shartemmpt at SSGI
// A LOT of help was given me by AlucardDH, creator of UberRT <3

#include "ReShade.fxh"
#include "soupcan_includes/random.fxh"
#include "soupcan_includes/pos.fxh"

// CONSTS
#define pi 3.14159265359
#define twoPi 6.28318530718
#define halfPi 1.57079632679

uniform int steps <> = 10;
uniform float thickness <> = 0.01;

struct Ray {
	float3 origin;
	float3 direction;
};

float3 at(Ray ray, float a) {
	return ray.origin + ray.direction * a;
}

float3 GetScreenSpaceNormal(float2 texcoord) {
    float3 offset = float3(BUFFER_PIXEL_SIZE, 0.0);
    float2 posCenter = texcoord.xy;
    float2 posNorth  = posCenter - offset.zy;
    float2 posEast   = posCenter + offset.xz;

    float3 vertCenter = float3(posCenter - 0.5, 1) * ReShade::GetLinearizedDepth(posCenter);
    float3 vertNorth  = float3(posNorth - 0.5,  1) * ReShade::GetLinearizedDepth(posNorth);
    float3 vertEast   = float3(posEast - 0.5,   1) * ReShade::GetLinearizedDepth(posEast);

    return normalize(cross(vertCenter - vertNorth, vertCenter - vertEast));
}

float4 trace(float3 normal, float3 pos, float3 direction, float2 texcoord) {
	Ray ray;
	ray.origin = pos;
	ray.direction = direction;
	
	for (int i = 1; i < steps; i++) {
		float3 step = at(ray, i);
		float3 screenPos = getScreenPosition(step);
		if (ReShade::GetLinearizedDepth(screenPos) < screenPos.z && screenPos.z < ReShade::GetLinearizedDepth(screenPos.z) + thickness) return tex2D(ReShade::BackBuffer, screenPos);
	}
	return 0;
}

float4 main(float4 vpos : SV_Position, float2 texcoord : TEXCOORD) : SV_Target {
	float depth = ReShade::GetLinearizedDepth(texcoord);
	uint ss = getPixelIndex(texcoord, BUFFER_SCREEN_SIZE);
	float3 randDir = randomVec3(texcoord, ss);
	float3 normal = GetScreenSpaceNormal(texcoord);
	if (dot(normal, randDir) < 0) randDir = -randDir;
	float3 pos = getWorldPosition(texcoord, depth);
	return trace(normal, pos, randDir, texcoord);
	//return float4(tex2Dfetch(ReShade::BackBuffer, vpos.xy));
}

technique SoupCanTraceRTGI {
	pass main {
		VertexShader = PostProcessVS;
		PixelShader = main;
	}
}