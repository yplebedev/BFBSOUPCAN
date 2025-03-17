/// My shartemmpt at SSGI
// Aknowledgements:
// A LOT of help was given me by AlucardDH, creator of UberRT <3
// I was set on the right path and consistently helped by Zenteon

#include "ReShade.fxh"
#include "soupcan_includes/random.fxh"
#include "soupcan_includes/pos.fxh"
#include "soupcan_includes/color_ops.fxh"

// CONSTS
#define pi 3.14159265359
#define twoPi 6.28318530718
#define halfPi 1.57079632679

float4 sampleBBlin(float2 texcoord) {
	return sRGBtoLinear(tex2D(ReShade::BackBuffer, texcoord));
}

uniform int steps <ui_type = "slider"; ui_min = 0.0; ui_max = 10.0;> = 10;
uniform float thickness <ui_type = "slider"; ui_min = 0.0; ui_max = 10.0;> = 0.01;
uniform float range <ui_type = "slider"; ui_min = 0.0; ui_max = 10.0;> = 1.0;
uniform float strength <ui_type = "slider"; ui_min = 0.0; ui_max = 10.0;> = 4.0;
uniform bool debug <> = false;

//uniform float3 lightPos <> = float3(0.0, 0.0, 0.0);

uniform int framecount < source = "framecount"; >;




float4 calculateLighting(float3 normal, float3 pos, float3 lightPos, float4 color) {
	float3 toLight = lightPos - pos;
	float dist = length(toLight);
	float light = max(dot(normalize(toLight), normal), 0.0) / (dist * dist + 1.0) * strength;
	float3 colored = color.rgb * light;
	return float4(colored, color.a);
}

float4 main(float4 vpos : SV_Position, float2 texcoord : TEXCOORD) : SV_Target {
	float4 input = sampleBBlin(texcoord);
	float depth = GetDepth(texcoord);
	if (1.0 - depth < 0.00001) return linearTosRGB(input);
	uint ss = getPixelIndex(texcoord, BUFFER_SCREEN_SIZE);
	uint ss2 = ss * 2;
	float3 randDir = randomVec3(texcoord, ss);
	float3 normal = GetWorldSpaceNormal(texcoord);
	float3 pos = getWorldPosition(texcoord, depth);
	
	float2 randF2 = float2(randomValue(ss), randomValue(ss2));
	float2 samplePos2D = float2(randF2.x * cos(randF2.y), randF2.y * cos(randF2.x));
	float3 samplePos3D = getWorldPosition(samplePos2D, GetDepth(samplePos2D + texcoord));
	
	float4 light = calculateLighting(-normal, pos, samplePos3D, sampleBBlin(samplePos2D + texcoord));
	if (!debug) {
		float4 x = linearTosRGB(saturate(input * light) * strength + input);
		return x / (x + 1);
	}
	return linearTosRGB(light + 0.5);
}

technique SoupCanTraceRTGI {
	pass main {
		VertexShader = PostProcessVS;
		PixelShader = main;
	}
}