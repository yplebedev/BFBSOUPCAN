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
#define zfw Zenteon 
#include "soupcan_includes/FrameworkResources.fxh"

#ifndef PI
	#define PI 3.14159265358979
#endif

#ifndef RENDER_MULT_X
	#define RENDER_MULT_X 0.5
#endif

#ifndef RENDER_MULT_Y
	#define RENDER_MULT_Y 0.5
#endif

texture bnt <source = "stbn.png";> {Width = 1024; Height = 1024; Format = R8; };
sampler bn { Texture = bnt; };


texture AOTex { Width = BUFFER_WIDTH * RENDER_MULT_X; Height = BUFFER_HEIGHT * RENDER_MULT_Y; Format = RGBA16F; };
sampler AOS { Texture = AOTex;
	MagFilter = POINT;
	MinFilter = POINT;
	MipFilter = POINT; };
texture AOTexPrev { Width = BUFFER_WIDTH * RENDER_MULT_X; Height = BUFFER_HEIGHT * RENDER_MULT_Y; Format = RGBA16F; };
sampler AOSPrev { Texture = AOTexPrev; 
	MagFilter = POINT;
	MinFilter = POINT;
	MipFilter = POINT; };
	
texture BBQ { Width = BUFFER_WIDTH / 4; Height = BUFFER_HEIGHT / 4; Format = RGBA16F; };
sampler BBQs { Texture = BBQ; 
	MagFilter = POINT;
	MinFilter = POINT;
	MipFilter = POINT; };

uniform int framecount < source = "framecount"; >;

uniform float historySize <ui_type = "slider"; ui_label = "Frame Blending"; ui_tooltip = "Affects the noise over update speed and ghosting ratios. This can be higher on higher FPS. 0 is no accumulation, and the closer to 1 the more previous results affect the image."; ui_min = 0.0; ui_max = 0.999;> = 0.8; 
uniform bool debug <ui_label = "Debug view";> = false;
uniform float strength <ui_type = "slider"; ui_label = "Strength"; ui_tooltip = "How much AO affects the input colors. Use conservativly."; ui_min = 0.0; ui_max = 4.0;> = 1.0; 

#define c_phi 1000.0 // Color avoiding
#define n_phi 0.3 // Normals
#define p_phi 1.0 // Depth

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

#define SLICES 6

#define STEPS 12

#define SECTORS 3

#define R 12000.0

#define R_MAX_CLAMP 13000.0

#define FAR_CLIP (RESHADE_DEPTH_LINEARIZATION_FAR_PLANE-1)

#define THICKNESS 2.0

// Modified from https://www.shadertoy.com/view/4djSRW
// It's called once per pixel.
float2 fastHash(float2 p) {
	float3 p3 = float3(p * BUFFER_SCREEN_SIZE, framecount & 0xFFFF);
	p3 = frac(p3 * float3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx + 33.33);
    return frac((p3.xy + p3.yz) * p3.zx);
}

float2 getTemporalOffset() {
	return float2(framecount % 8, (framecount >> 3) % 8);
}

// also vpos
float2 stbn(float2 p) {
	#define xyOffset float2(5, 7)
	return float2(tex2Dfetch(bn, (p % 64) + getTemporalOffset() * 64).x,
					tex2Dfetch(bn, ((p + xyOffset) % 64) + getTemporalOffset() * 64).x);
	// https://developer.nvidia.com/blog/rendering-in-real-time-with-spatiotemporal-blue-noise-textures-part-2/
	/* One way to do this is to read the texture at some fixed offset. For instance, if you read the first value at (pixelX, pixelY) % textureSize, you might read the second value at (pixelX+5, pixelY+7) % textureSize. This essentially gives you an uncorrelated spatiotemporal blue noise value, just as if you had a second spatiotemporal blue noise texture you were reading from.

The reason this works is because blue noise textures have correlation only over short distances. At long distances, the values are uncorrelated, as shown in Figure 7.

*/
}




float4 atrous(sampler input, float2 texcoord, float level) {
	float4 noisy = tex2D(input, texcoord);
	float3 normal = zfw::getNormal(texcoord);
	float depth = ReShade::GetLinearizedDepth(texcoord);
	if (1.0 - depth < 0.00001) return 1.0;
	float3 pos = zfw::uvToView(texcoord);
	
	float4 sum = 0.0;
	float2 step = ReShade::PixelSize;
	
	
	float cum_w = 0.0;
	[unroll]
	for (int i = 0; i < 25; i++) {
		float2 uv = texcoord + offset[i] * step * exp2(level);

		float4 ctmp = tex2Dlod(input, float4(uv, 0.0, 0.0));
		float4 t = noisy - ctmp;
		
		float dist2 = dot(t, t);
		float c_w = min(exp(-(dist2)/c_phi), 1.0);
		
		float3 ntmp = zfw::getNormal(uv);
		t = normal - ntmp;
		dist2 = max(dot(t, t), 0.0);
		float n_w = min(exp(-dist2 / n_phi), 1.0);
		
		float3 ptmp = zfw::uvToView(uv);
		t = pos - ptmp;
		dist2 = dot(t, t);
		float p_w = min(exp(-dist2 / p_phi), 1.0);
		
		float weight = c_w * n_w * p_w;
		sum += ctmp * weight * kernel[i];
		cum_w += weight * kernel[i];
	}
	return sum/cum_w;
}


uint sliceSteps(float3 positionVS, float3 V, float2 start, float2 rayDir, float t, float step, float samplingDirection, float N, uint bitfield, in out float3 gi, float3 n) {
    for (uint i = 0; i < STEPS; i++, t += step) {
        float2 samplePos = clamp(start + t * rayDir, 1, BUFFER_SCREEN_SIZE - 2);
        float3 di = tex2Dfetch(BBQs, samplePos / 4).rgb;
        float d = zfw::sampleDepth(samplePos.xy / BUFFER_SCREEN_SIZE , 0.0);
        float3 samplePosVS = zfw::uvzToView(float3(samplePos.xy / BUFFER_SCREEN_SIZE, d));
        float3 delta = samplePosVS - positionVS;
        float3 lj = normalize(delta);
        float falloff = rcp(dot(delta, delta) + exp2(-32));
	
	    float2 fb = acos(float2(dot(lj, V), dot(normalize(delta - V * THICKNESS), V)));
	    fb = saturate(((samplingDirection * -fb) - N + PI/2) / PI);
	    fb = samplingDirection >= 0 ? fb.yx : fb.xy;

   	 uint a = (fb.x * SECTORS);
    	uint b = ceil((fb.y - fb.x) * SECTORS);
    	uint prevBitfield = bitfield;
    	bitfield |= ((1 << b) - 1) << a;
    	float3 normal = zfw::sampleNormal(samplePos / BUFFER_SCREEN_SIZE, 0);
    	gi += countbits(bitfield & ~prevBitfield) 
		* di * saturate(dot(n, lj)) * saturate(dot(normal, -normalize(delta)))
		* dot(samplePosVS, samplePosVS) * (t / float(STEPS))
		* falloff / SECTORS;
    }
    return bitfield;
}

float4 saveBB(float4 vpos : SV_Position, float2 uv : TEXCOORD) : SV_Target {
	float3 bb = zfw::toneMapInverse(tex2Dfetch(ReShade::BackBuffer, vpos.xy * 4).rgb, 20.0);
	float3 cached = tex2Dfetch(AOS, vpos.xy / float2(RENDER_MULT_X, RENDER_MULT_Y)).rgb * strength * zfw::getAlbedo(uv).rgb;
	return float4(bb + cached, 1.0);
}

float4 gtao(float2 uv, float2 vpos) {
	float2 random = stbn(vpos);
	
	float2 start = uv * BUFFER_SCREEN_SIZE;
	float3 positionVS = zfw::uvToView(uv);
	positionVS.z *= 0.9999; // Move center pixel towards camera a bit.

	float ao = 0.0;
	float3 gi = 0.0;
	
	float3 V = normalize(-positionVS);
	float3 normalVS = zfw::getNormal(uv);

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
		//float offset = step;
		aoBF = sliceSteps(positionVS, V, start, direction, offset, step, 1, N, aoBF, gi, normalVS);
		aoBF = sliceSteps(positionVS, V, start, -direction, offset, step, -1, N, aoBF, gi, normalVS);

		ao += float(countbits(aoBF));
	}
	ao = 1.0 - ao / (float(SECTORS) * SLICES);
	return float4(gi / (STEPS * SLICES), positionVS.z > FAR_CLIP || ao < -0.001 ? 1.0 : ao);
}

float4 main(float4 vpos : SV_Position, float2 uv : TEXCOORD) : SV_Target {
	float3 MV = zfw::getVelocity(uv);
	return lerp(gtao(uv, vpos.xy), tex2D(AOSPrev, uv + MV.xy), historySize * MV.z);
}


float4 display(float4 vpos : SV_Position, float2 uv : TEXCOORD) : SV_Target {
	float4 denoised = atrous(AOS, uv, 0);
	float3 bb = zfw::getBackBuffer(uv);
	float3 lighting = denoised.rgb * denoised.a * zfw::getAlbedo(uv);
	if (!debug)
		return float4(zfw::toneMap(zfw::toneMapInverse(bb, 20.0) + lighting, 20.0), 1.0);
	else
		return float4(zfw::toneMap((denoised.rgb + (0.05 * denoised.a)) * strength, 20.0), 1.0);
}

float4 cache(float4 vpos : SV_Position, float2 uv : TEXCOORD) : SV_Target {
	return tex2Dfetch(AOS, vpos.xy);
}

technique HBIL {
	pass BB {
		VertexShader = PostProcessVS;
		PixelShader = saveBB;
		RenderTarget = BBQ;
	}
	pass Main { 
		VertexShader = PostProcessVS;
		PixelShader = main;
		RenderTarget = AOTex;
	}
	pass Display {
		VertexShader = PostProcessVS;
		PixelShader = display;
	}
	pass Cache {
		VertexShader = PostProcessVS;
		PixelShader = cache;
		RenderTarget = AOTexPrev;
	}
}
