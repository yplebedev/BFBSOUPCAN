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


texture AOTex { Width = BUFFER_WIDTH * RENDER_MULT_X; Height = BUFFER_HEIGHT * RENDER_MULT_Y; Format = R8; };
sampler AOS { Texture = AOTex;
	MagFilter = POINT;
	MinFilter = POINT;
	MipFilter = POINT; };
texture AOTexPrev { Width = BUFFER_WIDTH * RENDER_MULT_X; Height = BUFFER_HEIGHT * RENDER_MULT_Y; Format = R8; };
sampler AOSPrev { Texture = AOTexPrev; 
	MagFilter = POINT;
	MinFilter = POINT;
	MipFilter = POINT; };

uniform int framecount < source = "framecount"; >;

uniform float historySize <ui_type = "slider"; ui_label = "Frame Blending"; ui_tooltip = "Affects the noise over update speed and ghosting ratios. This can be higher on higher FPS. 0 is no accumulation, and the closer to 1 the more previous results affect the image."; ui_min = 0.0; ui_max = 0.999;> = 0.8; 
uniform bool debug <ui_label = "Debug view";> = false;
uniform float strength <ui_type = "slider"; ui_label = "Strength"; ui_tooltip = "How much AO affects the input colors. Use conservativly."; ui_min = 0.0; ui_max = 1.0;> = 1.0; 
uniform float power <ui_type = "slider"; ui_label = "Sampling bias power"; ui_tooltip = "A bias to sample closer-by stuff more. Increasing increases effective radius and occlusion from small objects."; ui_min = 1.0; ui_max = 3.0;> = 1.0; 
uniform float R <ui_type = "slider"; ui_label = "Radius"; ui_tooltip = "Increases the effect scale, but lowers effective quality and lowers cache coherency."; ui_min = 50.0; ui_max = 50000.0;> = 3000.0; 
uniform float THICKNESS <ui_type = "slider"; ui_label = "Thickness"; ui_tooltip = "SCVBAO uses... you guessed it, VBAO. It allows for a thickness heuristic. Don't set this too high or low!"; ui_min = 0.0; ui_max = 4.0;> = 2.0; 

#define c_phi 1 // Color avoiding
#define n_phi 1 // Normals
#define p_phi 1 // Depth

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

#ifndef SCVBAO_SLICES
	#define SCVBAO_SLICES 2
#endif

#ifndef SCVBAO_STEPS
		#define SCVBAO_STEPS 7
#endif

#define SECTORS 32

//#define R 3000.0

#define R_MAX_CLAMP 5000

#define FAR_CLIP (RESHADE_DEPTH_LINEARIZATION_FAR_PLANE-1)

//#define THICKNESS 2

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
	
}

uint sliceSteps(float3 positionVS, float3 V, float2 start, float2 rayDir, float t, float step, float samplingDirection, float N, uint bitfield) {
    for (uint i = 0; i < SCVBAO_STEPS; i++, t += step) {
        float2 samplePos = clamp(start + t * rayDir, 1, BUFFER_SCREEN_SIZE - 2);
        samplePos = floor(samplePos) + 0.5;
        float3 samplePosVS = zfw::uvToView(samplePos.xy / BUFFER_SCREEN_SIZE);
        float3 delta = samplePosVS - positionVS;
	
	    float2 fb = acos(float2(dot(normalize(delta), V), dot(normalize(delta + THICKNESS * normalize(samplePosVS)), V)));
	    fb = saturate(((samplingDirection * -fb) - N + PI/2) / PI);
	    fb = fb.x > fb.y ? fb.yx : fb;
	    fb = smoothstep(0., 1., fb); // cosine lobe for AO. look: cdf.
	    
   	 uint a = round(fb.x * SECTORS);
    	uint b = round((fb.y - fb.x) * SECTORS);
    	bitfield |= ((1 << b) - 1) << a;
    	t *= power;
    }
    return bitfield;
}


float gtao(float2 uv, float2 vpos) {
	float2 random = stbn(vpos);
	
	float2 start = vpos;
	float3 positionVS = zfw::uvToView(uv);
	positionVS.z *= 0.9999; // Move center pixel towards camera a bit.

	float ao = 0.0;
	
	float3 V = normalize(-positionVS);
	float3 normalVS = zfw::getNormal(uv);

    float step = max(1.0, clamp(R / positionVS.z, SCVBAO_STEPS, R_MAX_CLAMP) / (SCVBAO_STEPS + 1.0));
		
	for(float slice = 0.0; slice < 1.0; slice += 1.0 / SCVBAO_SLICES) {
		float phi = PI * frac(slice + random.x);
		float2 direction = float2(cos(phi), sin(phi));
		
		float3 directionF3 = float3(direction, 0.0);
		float3 oDirV = directionF3 - dot(directionF3, V) * V;
		float3 sliceN = cross(directionF3, V);
		float3 projN = normalVS - sliceN * dot(normalVS, sliceN);
		float cosN = saturate(dot(projN, V) / length(projN));
		float signN = -sign(dot(projN, oDirV));

		float N = signN * acos(cosN);
		
		uint aoBF = 0;
		float offset = max(random.y * step, length(BUFFER_PIXEL_SIZE));
		aoBF = sliceSteps(positionVS, V, start, direction, offset, step, 1, N, aoBF);
		aoBF = sliceSteps(positionVS, V, start, -direction, offset, step, -1, N, aoBF);

		ao += float(countbits(aoBF));
	}
	ao = 1.0 - ao / (float(SECTORS) * SCVBAO_SLICES);
	return positionVS.z > FAR_CLIP || ao < -0.001 ? 1.0 : ao;
}

float4 main(float4 vpos : SV_Position, float2 uv : TEXCOORD) : SV_Target {
	const float ao = pow(gtao(uv, vpos.xy), 1);
	return float4(ao, ao, ao, 1.0);
}


technique SCAO {
	pass Main { 
		VertexShader = PostProcessVS;
		PixelShader = main;
	}
}
