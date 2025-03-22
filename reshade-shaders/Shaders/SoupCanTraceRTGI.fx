/// My shartemmpt at SSGI
// Aknowledgements:
// A LOT of help was given me by AlucardDH, creator of UberRT <3
// I was set on the right path and consistently helped by Zenteon
// The quasirandom sequence is an extention of R2 by Marty, creator of the famous RTGI shader.

#include "ReShade.fxh"
#include "soupcan_includes/random.fxh"
#include "soupcan_includes/pos.fxh"
#include "soupcan_includes/color_ops.fxh"

#define BUFFER_SIZE int2(BUFFER_WIDTH,BUFFER_HEIGHT)
#define NOISE_TEX_NAME "QuarkBN_SOUPGI.png"

texture bn < source = NOISE_TEX_NAME; > { Width = 512; Height = 512; Format = RGBA8; };
sampler sbn { Texture = bn; };


// CONSTS
#define pi 3.14159265359
#define twoPi 6.28318530718
#define halfPi 1.57079632679


#ifndef RENDER_SCALE
	#define RENDER_SCALE 0.5
	#define BUFFER_SCALED_SIZE = BUFFER_SIZE * RENDER_SCALE
#endif

#ifndef SOUPCANGI_SPP
	#define SOUPCANGI_SPP 1
#endif

texture DN0 { Width = BUFFER_WIDTH * RENDER_SCALE; Height = BUFFER_HEIGHT * RENDER_SCALE; Format = RGBA16F; };
texture DN1 { Width = BUFFER_WIDTH * RENDER_SCALE; Height = BUFFER_HEIGHT * RENDER_SCALE; Format = RGBA16F; };

sampler SDN0 { Texture = DN0; };
sampler SDN1 { Texture = DN1; };

float4 sampleBBlin(float2 texcoord) {
	return sRGBtoLinear(tex2D(ReShade::BackBuffer, texcoord));
}

uniform float strength <ui_type = "slider"; ui_min = 0.0; ui_max = 10.0;> = 4.0;
uniform bool debug <> = false;

//uniform float3 lightPos <> = float3(0.0, 0.0, 0.0);

uniform int framecount < source = "framecount"; >;

texture intermediateLight { Width = BUFFER_WIDTH * RENDER_SCALE; Height = BUFFER_HEIGHT * RENDER_SCALE; Format = RGBA16F; };
sampler intermediateLightSampler { Texture = intermediateLight; };

uniform float c_phi <ui_type = "slider"; ui_min = 0.1; ui_max = 10.0;> = 0.1;
uniform float n_phi <ui_type = "slider"; ui_min = 0.1; ui_max = 10.0;> = 0.275;
uniform float p_phi <ui_type = "slider"; ui_min = 0.1; ui_max = 10.0;> = 0.161;

uniform float kernel[25] <ui_title = "GAUSS WHEIGHTS - DO NOT EDIT";> = {
    1.0/256.0, 1.0/64.0,  3.0/128.0, 1.0/64.0,  1.0/256.0,
    1.0/64.0,  1.0/16.0,  3.0/32.0,  1.0/16.0,  1.0/64.0,
    3.0/128.0, 3.0/32.0,  9.0/64.0,  3.0/32.0,  3.0/128.0,
    1.0/64.0,  1.0/16.0,  3.0/32.0,  1.0/16.0,  1.0/64.0,
    1.0/256.0, 1.0/64.0,  3.0/128.0, 1.0/64.0,  1.0/256.0
};
uniform float2 offset[25] <ui_title = "OFFSETS - DO NOT EDIT";> = {
    float2(-2.0, -2.0), float2(-1.0, -2.0), float2(0.0, -2.0), float2(1.0, -2.0), float2(2.0, -2.0),
    float2(-2.0, -1.0), float2(-1.0, -1.0), float2(0.0, -1.0), float2(1.0, -1.0), float2(2.0, -1.0),
    float2(-2.0,  0.0), float2(-1.0,  0.0), float2(0.0,  0.0), float2(1.0,  0.0), float2(2.0,  0.0),
    float2(-2.0,  1.0), float2(-1.0,  1.0), float2(0.0,  1.0), float2(1.0,  1.0), float2(2.0,  1.0),
    float2(-2.0,  2.0), float2(-1.0,  2.0), float2(0.0,  2.0), float2(1.0,  2.0), float2(2.0,  2.0)
};


// ToDo: ???? major refactor in order, this is about as clean as ZN's first repo lmao.
float4 calculateLighting(float3 normal, float3 pos, float3 lightPos, float4 color, float3 normalAtSampled, float2 texcoord) {
	float3 toLight = lightPos - pos;
	float3 fromLight = normalize(-toLight);
	float dist = length(toLight);
	float light = max(dot(normalize(toLight), normal), 0.0) / (dist * dist * 0.01 + 1.0) * strength;
	//light /= dot(pos, pos);
	light *= dot(lightPos, lightPos);
	float weight = saturate(dot(-normalAtSampled, fromLight));
	float3 colored = color.rgb * light * weight;
	return float4(colored, color.a);
}

float4 main(float4 vpos : SV_Position, float2 texcoord : TEXCOORD) : SV_Target {
	float4 input = inverseTonemapLottes(sampleBBlin(texcoord));
	float depth = GetDepth(texcoord);
	if (1.0 - depth < 0.00001) return linearTosRGB(input);
	float3 normal = -GetWorldSpaceNormal(texcoord);
	float3 pos = getWorldPosition(texcoord, depth);
	float2 seed = tex2Dfetch(sbn, vpos.xy % 512).rg;
	
	float4 light = 0;
	[unroll]
	for (int i = 0; i < SOUPCANGI_SPP; i++) {
		float2 randF2 = r2_modified(i, seed);
		float3 normalAtSampled = GetWorldSpaceNormal(randF2);
		float3 samplePos3D = getWorldPosition(randF2, GetDepth(randF2));
		
		light += calculateLighting(normal, pos, samplePos3D, inverseTonemapLottes(sampleBBlin(randF2)), normalAtSampled, texcoord);
	}
	return light / SOUPCANGI_SPP;
}

float4 upscale(float4 vpos, float2 texcoord, sampler s, float level) {
// https://jo.dreggn.org/home/2010_atrous.pdf
	float4 noisy = tex2D(s, texcoord);
	float3 normal = GetWorldSpaceNormal(texcoord);
	float depth = GetDepth(texcoord);
	if (1.0 - depth < 0.00001) return 0;
	float3 pos = getWorldPosition(texcoord, depth);
	
	float4 sum = 0.0;
	float2 step = ReShade::PixelSize;
	
	
	float cum_w = 0.0;
	[unroll]
	for (int i = 0; i < 25; i++) {
		float2 uv = texcoord + offset[i] * step * exp2(level);
		
		float3 ctmp = tex2Dlod(s, float4(uv, 0.0, 0.0)).rgb;
		float3 t = noisy.rgb - ctmp;
		
		float dist2 = dot(t, t);
		float c_w = min(exp(-(dist2)/c_phi), 1.0);
		
		float3 ntmp = GetWorldSpaceNormal(uv);
		t = normal - ntmp;
		dist2 = max(dot(t, t), 0.0);
		float n_w = min(exp(-dist2 / n_phi), 1.0);
		
		float depth = GetDepth(uv);
		float3 ptmp = getWorldPosition(uv, depth);
		t = pos - ptmp;
		dist2 = dot(t, t);
		float p_w = min(exp(-dist2 / p_phi), 1.0);
		
		float weight = c_w * n_w * p_w;
		sum += ctmp * weight * kernel[i];
		cum_w += weight * kernel[i];
	}
	
	float4 light = sum/cum_w;
	return light;
}


float4 upscalePass0(float4 vpos : SV_Position, float2 texcoord : TEXCOORD) : SV_Target {
	return upscale(vpos, texcoord, intermediateLightSampler, 0);
}

float4 upscalePass1(float4 vpos : SV_Position, float2 texcoord : TEXCOORD) : SV_Target {
	return upscale(vpos, texcoord, SDN1, 1);
}

float4 upscalePass2(float4 vpos : SV_Position, float2 texcoord : TEXCOORD) : SV_Target {
	return upscale(vpos, texcoord, SDN0, 2);
}

float4 blend(float4 vpos : SV_Position, float2 texcoord : TEXCOORD) : SV_Target {
	return tonemapLottes(tex2Dfetch(SDN1, vpos.xy * RENDER_SCALE));
}

technique SoupCanTraceRTGI {
	pass light {
		VertexShader = PostProcessVS;
		PixelShader = main;
		RenderTarget = intermediateLight;
	}
	pass upscale0 {
		VertexShader = PostProcessVS;
		PixelShader = upscalePass0;
		RenderTarget = DN1;
	}
	pass upscale1 {
		VertexShader = PostProcessVS;
		PixelShader = upscalePass1;
		RenderTarget = DN0;
	}
	pass upscale2 {
		VertexShader = PostProcessVS;
		PixelShader = upscalePass2;
		RenderTarget = DN1;
	}
	pass blend {
		VertexShader = PostProcessVS;
		PixelShader = blend;
	}
}