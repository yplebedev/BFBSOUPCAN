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


texture AO1 { Width = BUFFER_WIDTH / 2; Height = BUFFER_HEIGHT / 2; Format = R8; };
sampler sAO1 { Texture = AO1;
	MagFilter = POINT;
	MinFilter = POINT;
	MipFilter = POINT; };
texture AO2 { Width = BUFFER_WIDTH / 2; Height = BUFFER_HEIGHT / 2; Format = R8; };
sampler sAO2 { Texture = AO2; 
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
	#define SCVBAO_SLICES 1
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

float4 atrous(sampler input, float2 texcoord, float level) {
	float4 noisy = tex2D(input, texcoord);
	float3 normal = zfw::getNormal(texcoord);
	float depth = ReShade::GetLinearizedDepth(texcoord);
	if (1.0 - depth < 0.00001) return 1.0;
	float3 pos = zfw::uvToView(texcoord);
	
	float4 sum = 0.0;
	float2 step = ReShade::PixelSize * 2;
	
	
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

sampler lowN { Texture = zfw::tLowNormal; };
sampler highN{ Texture = zfw::tNormal; };
// From papadanku.github.io
// https://www.semanticscholar.org/paper/Multistep-joint-bilateral-depth-upsampling-Riemens-Gangwal/1ddf9ad017faf63b04778c1ddfc2330d64445da8
float4 JointBilateralUpsample(
   sampler Image, // This should be 1/2 the size as GuideHigh
   sampler GuideLow, // This should be 1/2 the size as GuideHigh
   sampler GuideHigh, // This should be 2/1 the size as Image and GuideLow
   float2 Tex ) {
	   // Initialize variables
	   float2 PixelSize = ldexp(fwidth(Tex.xy), 1.0) * 0.5;
	   float4 GuideHighSample = tex2D(GuideHigh, Tex);
	   float4 BilateralSum = 0.0;
	   float4 WeightSum = 0.0;
	
	   [unroll]
	   for (int dx = -1; dx <= 1; ++dx) {
	      [unroll]
	      for (int dy = -1; dy <= 1; ++dy) {
	         // Calculate offset
	         float2 Offset = float2(float(dx), float(dy));
	         float2 OffsetTex = Tex + (Offset * PixelSize);
	
	         // Sample image and guide
	         float4 ImageSample = tex2Dlod(Image, float4(OffsetTex, 0.0, 0.0));
	         float4 GuideLowSample = tex2D(GuideLow, OffsetTex);
	
	         // Calculate weight
	         float3 Delta = GuideLowSample.xyz - GuideHighSample.xyz;
	         float DotDD = dot(Delta, Delta);
	         float Weight = (DotDD > 0.0) ? 1.0 / DotDD : 1.0;
	
	         BilateralSum += (ImageSample * Weight);
	         WeightSum += Weight;
      }
   }

   return BilateralSum / WeightSum;
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
	
	float2 start = vpos * 2.;
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
	const float ao = pow(gtao(uv, vpos.xy / 2), 1);
	return float4(ao, ao, ao, 1.0);
}

float4 denoise1(float4 vpos : SV_Position, float2 uv : TEXCOORD) : SV_Target {
	return atrous(sAO1, uv, 1).xxxx;
}

float4 upscale(float4 vpos : SV_Position, float2 uv : TEXCOORD) : SV_Target {
	return JointBilateralUpsample(sAO2, lowN, highN, uv);
}

technique SCAO {
	pass Main { 
		VertexShader = PostProcessVS;
		PixelShader = main;
		RenderTarget = AO1;
	}
	pass Denoise {
		VertexShader = PostProcessVS;
		PixelShader = denoise1;
		RenderTarget = AO2;
	}
	pass Upscale {
		VertexShader = PostProcessVS;
		PixelShader = main;
	}
}
