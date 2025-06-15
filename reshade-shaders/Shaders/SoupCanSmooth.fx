// Cross bilateral filter approximation using the A-Trous filter.
// Sky is black. But we skip it anyway.
#include "ReShade.fxh"
#include "soupcan_includes/pos.fxh"

#define c_phi 1 // Color avoiding
#define n_phi 1 // Normals
#define p_phi 1 // Depth

texture DN0 { Width = BUFFER_WIDTH; Height = BUFFER_HEIGHT; Format = RGBA8; };
texture DN1 { Width = BUFFER_WIDTH; Height = BUFFER_HEIGHT; Format = RGBA8; };

sampler SDN0 { Texture = DN0; };
sampler SDN1 { Texture = DN1; };


uniform float kernel[25] <hidden = true;> = {
    1.0/256.0, 1.0/64.0,  3.0/128.0, 1.0/64.0,  1.0/256.0,
    1.0/64.0,  1.0/16.0,  3.0/32.0,  1.0/16.0,  1.0/64.0,
    3.0/128.0, 3.0/32.0,  9.0/64.0,  3.0/32.0,  3.0/128.0,
    1.0/64.0,  1.0/16.0,  3.0/32.0,  1.0/16.0,  1.0/64.0,
    1.0/256.0, 1.0/64.0,  3.0/128.0, 1.0/64.0,  1.0/256.0
};
uniform float2 offset[25] <hidden = true;> = {
    float2(-2.0, -2.0), float2(-1.0, -2.0), float2(0.0, -2.0), float2(1.0, -2.0), float2(2.0, -2.0),
    float2(-2.0, -1.0), float2(-1.0, -1.0), float2(0.0, -1.0), float2(1.0, -1.0), float2(2.0, -1.0),
    float2(-2.0,  0.0), float2(-1.0,  0.0), float2(0.0,  0.0), float2(1.0,  0.0), float2(2.0,  0.0),
    float2(-2.0,  1.0), float2(-1.0,  1.0), float2(0.0,  1.0), float2(1.0,  1.0), float2(2.0,  1.0),
    float2(-2.0,  2.0), float2(-1.0,  2.0), float2(0.0,  2.0), float2(1.0,  2.0), float2(2.0,  2.0)
};

//				BackBuffer        uv              iteration
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
	for (int i = 0; i < 25; i++) {
		float2 uv = texcoord + offset[i] * step * exp2(level);

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

float4 d0(float4 vpos : SV_Position, float2 texcoord : TEXCOORD) : SV_Target {
	return float4(atrous(ReShade::BackBuffer, texcoord, 0), 1.0);
}

float4 d1(float4 vpos : SV_Position, float2 texcoord : TEXCOORD) : SV_Target {
	return float4(atrous(SDN0, texcoord, 1), 1.0);
}

float4 d2(float4 vpos : SV_Position, float2 texcoord : TEXCOORD) : SV_Target {
	return float4(atrous(SDN1, texcoord, 2), 1.0);
}

technique SoupCanSmooth {
	pass {
		PixelShader = d0;
		VertexShader = PostProcessVS;
		RenderTarget = DN0;
	}
	pass {
		PixelShader = d1;
		VertexShader = PostProcessVS;
		RenderTarget = DN1;
	}
	pass {
		PixelShader = d2;
		VertexShader = PostProcessVS;
	}
}
