#include "ReShade.fxh"
#include "soupcan_includes/ao.fxh"

#ifndef DEBUG_SCAO
	#define DEBUG_SCAO false
#endif

#ifndef SCAO_RENDER_SCALE
	#define SCAO_RENDER_SCALE 0.5
	#define RCP_RENDER_SCALE rcp(SCAO_RENDER_SCALE)
#endif

#define c_phi 100
#define n_phi 1
#define p_phi 0.4


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



texture AORT { Width = BUFFER_WIDTH * SCAO_RENDER_SCALE; Height = BUFFER_HEIGHT * SCAO_RENDER_SCALE; Format = RGBA8; };
sampler sAORT { Texture = AORT; };

// some day ill prolly do something useful here.
float4 computeNormals(float4 vpos : SV_Position, float2 texcoord : TEXCOORD) : SV_Target {
	return getPackedWorldSpaceNormal(texcoord);
}

float3 atrous(sampler input, float2 texcoord, float level) {
	float4 noisy = tex2D(input, texcoord);
	float3 normal = getWorldSpaceNormal(texcoord);
	float depth = GetDepth(texcoord);
	if (1.0 - depth < 0.00001) return 0;
	float3 pos = getWorldPosition(texcoord, depth);
	
	float3 sum = 0.0;
	float2 step = ReShade::PixelSize;
	
	
	float cum_w = 0.0;
	[unroll]
	for (int i = 0; i < 25; i++) {// 										this might suck
		float2 uv = texcoord + offset[i] * step * exp2(level) * SCAO_RENDER_SCALE;
		//																   it is smoother :)
		float3 ctmp = tex2Dlod(input, float4(uv, 0.0, 0.0)).rgb;
		float3 t = noisy.rgb - ctmp;
		
		float dist2 = dot(t, t);
		float c_w = min(exp(-(dist2)/c_phi), 1.0);
		
		float3 ntmp = getWorldSpaceNormal(uv);
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
	return sum/cum_w;
}

float occlude(float4 vpos : SV_Position, float2 texcoord : TEXCOORD) : SV_Target {
	return gtao(texcoord, vpos.xy);
}

float4 blend(float4 vpos : SV_Position, float2 texcoord : TEXCOORD) : SV_Target {
	float4 input = tex2Dfetch(ReShade::BackBuffer, vpos.xy);
	float ao = atrous(sAORT, texcoord, 2).r;
	if (ao == 0) return input; // what the fuck.
	return float4(lerp(ao, input.rgb * ao, 1 - DEBUG_SCAO), 1.0);
}

technique SoupCanAO {
	pass AO {
		VertexShader = PostProcessVS;
		PixelShader = occlude;
		RenderTarget = AORT;
	}
	pass blend {
		VertexShader = PostProcessVS;
		PixelShader = blend;
	}
}