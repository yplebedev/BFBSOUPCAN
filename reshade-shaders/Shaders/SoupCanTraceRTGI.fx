/// My shartemmpt at SSGI
// A LOT of help was given me by AlucardDH, creator of UberRT <3

#include "ReShade.fxh"
#include "soupcan_includes/random.fxh"
#include "soupcan_includes/pos.fxh"
#include "soupcan_includes/color_ops.fxh"

// CONSTS
#define pi 3.14159265359
#define twoPi 6.28318530718
#define halfPi 1.57079632679
#define GetDepth ReShade::GetLinearizedDepth

float4 sampleBBlin(float2 texcoord) {
	return sRGBtoLinear(tex2D(ReShade::BackBuffer, texcoord));
}

uniform int steps <ui_type = "slider"; ui_min = 0.0; ui_max = 10.0;> = 10;
uniform float thickness <ui_type = "slider"; ui_min = 0.0; ui_max = 10.0;> = 0.01;
uniform float range <ui_type = "slider"; ui_min = 0.0; ui_max = 10.0;> = 1.0;
uniform float strength <ui_type = "slider"; ui_min = 0.0; ui_max = 10.0;> = 4.0;
uniform bool debug <> = false;

uniform int framecount < source = "framecount"; >;

namespace Rays {
	
	struct Ray {
		float3 origin;
		float3 direction;
	};
	
	float3 at(Ray ray, float a) {
		return ray.origin + ray.direction * a;
	}
	
	Ray offsetByNormal(Ray ray, float3 normal, float offset) {
		Ray rayOut;
		rayOut.origin = ray.origin + normal * offset;
		rayOut.direction = ray.direction + normal * offset; // Is this required?
		return rayOut;
	}
}

// ToDo: Refactor
namespace Hit {
	struct Hit {
		bool hit;
		float3 worldSpacePos;
		float3 worldSpaceNormal;
		float3 screenPos;
		float depth;
	};
	
	struct Interval {
		float3 start;
		float3 end;
	};
	
	float3 GetWorldSpaceNormal(float2 texcoord) {
	    float3 offset = float3(BUFFER_PIXEL_SIZE, 0.0);
	    float2 posCenter = texcoord.xy;
	    float2 posNorth  = posCenter - offset.zy;
	    float2 posEast   = posCenter + offset.xz;
		float3 vertCenter = getWorldPosition(posCenter,       GetDepth(posCenter));
	    float3 vertNorth  = getWorldPosition(posNorth,  GetDepth(posNorth));
	    float3 vertEast   = getWorldPosition(posEast,   GetDepth(posEast));
	    
	    return normalize(cross(vertCenter - vertNorth, vertCenter - vertEast));
	}

	
	Hit consider(Rays::Ray ray, float length) {
		Interval toConsider;
		toConsider.start = ray.origin;
		toConsider.end = Rays::at(ray, length);
		
		Hit hit;
		for (int i = 0; i < steps; i++) {
			hit.screenPos = getScreenPosition(toConsider.end);
			hit.depth = GetDepth(hit.screenPos.xy);
			hit.worldSpaceNormal = GetWorldSpaceNormal(hit.screenPos.xy);
			if (hit.depth < hit.screenPos.z && hit.depth + (thickness * hit.screenPos.z)> hit.screenPos.z) {
				toConsider.start = 
			} else {
				toConsider.
			}
		}
		
		return hit;
	}
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

float3 GetWorldSpaceNormal(float2 texcoord) {
    float3 offset = float3(BUFFER_PIXEL_SIZE, 0.0);
    float2 posCenter = texcoord.xy;
    float2 posNorth  = posCenter - offset.zy;
    float2 posEast   = posCenter + offset.xz;
	float3 vertCenter = getWorldPosition(posCenter,       GetDepth(posCenter));
    float3 vertNorth  = getWorldPosition(posNorth,  GetDepth(posNorth));
    float3 vertEast   = getWorldPosition(posEast,   GetDepth(posEast));
    
    return normalize(cross(vertCenter - vertNorth, vertCenter - vertEast));
}



float4 trace(float3 normal, float3 pos, float3 direction, float2 texcoord, uint seed) {
	Rays::Ray ray;
	ray.origin = pos;
	ray.direction = direction + normal;
	ray = Rays::offsetByNormal(ray, normal, 0.001);
	for (int i = 1; i < steps; i++) {
		seed = seed + framecount;
		float3 step = Rays::at(ray, i * (range + randomValue(seed)));
		float3 screenPos = getScreenPosition(step);
		float depth = GetDepth(screenPos.xy);
		float3 normalAtHit = GetWorldSpaceNormal(screenPos.xy);
		if (dot(normalAtHit, normal) < 0.0) return 0;
		if (depth < screenPos.z && depth + (thickness * screenPos.z)> screenPos.z) return inverseTonemapLottes(sampleBBlin(screenPos.xy));
	}
	return 0;
}

float4 main(float4 vpos : SV_Position, float2 texcoord : TEXCOORD) : SV_Target {
	float depth = ReShade::GetLinearizedDepth(texcoord);
	float4 input = sampleBBlin(texcoord);
	if (1.0 - depth < 0.00001) return linearTosRGB(input);
	uint ss = getPixelIndex(texcoord, BUFFER_SCREEN_SIZE);
	float3 randDir = randomVec3(texcoord, ss);
	float3 normal = GetWorldSpaceNormal(texcoord);
	//if (dot(-normal, randDir) < 0) randDir = -randDir;
	float3 pos = getWorldPosition(texcoord, depth);
	float4 traceRes = trace(-normal, pos, randDir, texcoord, ss);
	if (!debug) return linearTosRGB(saturate(input * traceRes) * strength + input);
	return linearTosRGB(traceRes + 0.5);
}

technique SoupCanTraceRTGI {
	pass main {
		VertexShader = PostProcessVS;
		PixelShader = main;
	}
}