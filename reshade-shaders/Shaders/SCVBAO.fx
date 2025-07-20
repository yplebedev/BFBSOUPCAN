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

uniform int framecount < source = "framecount"; >;


texture bnt <source = "stbn.png";> {Width = 1024; Height = 1024; Format = R8; };
sampler bn { Texture = bnt; };

/*
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
	
texture DN0 { Width = BUFFER_WIDTH * 0.5; Height = BUFFER_HEIGHT * 0.5; Format = R8; };
sampler sDN0 { Texture = DN0;
	MagFilter = POINT;
	MinFilter = POINT;
	MipFilter = POINT; };
texture DN1 { Width = BUFFER_WIDTH * 0.5; Height = BUFFER_HEIGHT * 0.5; Format = R8; };
sampler sDN1 { Texture = DN1; 
	MagFilter = POINT;
	MinFilter = POINT;
	MipFilter = POINT; };
	
texture Accum {  Width = BUFFER_WIDTH * RENDER_MULT_X; Height = BUFFER_HEIGHT * RENDER_MULT_Y; Format = R32U; };
sampler2D<uint> sAccum { Texture = Accum; };
texture AccumW {  Width = BUFFER_WIDTH * RENDER_MULT_X; Height = BUFFER_HEIGHT * RENDER_MULT_Y; Format = R32U; };
sampler2D<uint> sAccumW { Texture = AccumW; };

texture PrevN {  Width = BUFFER_WIDTH * RENDER_MULT_X; Height = BUFFER_HEIGHT * RENDER_MULT_Y; Format = RGBA8; };
sampler sPrevN { Texture = PrevN; };

uniform float historySize <ui_type = "slider"; ui_label = "Frame Blending"; ui_tooltip = "Affects the noise over update speed and ghosting ratios. This can be higher on higher FPS. 0 is no accumulation, and the closer to 1 the more previous results affect the image."; ui_min = 0.0; ui_max = 0.999;> = 0.8; 
uniform bool debug <ui_label = "Debug view";> = false;
uniform float strength <ui_type = "slider"; ui_label = "Strength"; ui_tooltip = "How much AO affects the input colors. Use conservativly."; ui_min = 0.0; ui_max = 1.0;> = 1.0; 

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
*/
#define SLICES 4

#define STEPS 10

#define SECTORS 32

#define R 50.0

#define R_MAX_CLAMP 3000

#define FAR_CLIP (RESHADE_DEPTH_LINEARIZATION_FAR_PLANE-1)

#define THICKNESS 0.4


// By Marty: https://www.martysmods.com/
float2 r2_modified(in int idx, in float2 seed) {
	return frac(seed + float(idx) * float2(0.245122333753, 0.430159709002));
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
	//. One way to do this is to read the texture at some fixed offset. For instance, if you read the first value at (pixelX, pixelY) % textureSize, you might read the second value at (pixelX+5, pixelY+7) % textureSize. This essentially gives you an uncorrelated spatiotemporal blue noise value, just as if you had a second spatiotemporal blue noise texture you were reading from.

//The reason this works is because blue noise textures have correlation only over short distances. At long distances, the values are uncorrelated, as shown in Figure 7.


}
/*
float atrous(sampler input, float2 texcoord, float level) {
	float4 noisy = tex2D(input, texcoord);
	float3 normal = zfw::getNormal(texcoord);
	float depth = ReShade::GetLinearizedDepth(texcoord);
	if (1.0 - depth < 0.00001) return 1.0;
	float3 pos = zfw::uvToView(texcoord);
	
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
	return sum.x/cum_w.x;
}

sampler lowN { Texture = zfw::tLowNormal; };
sampler highN{ Texture = zfw::tNormal; }; //ugly ahh nonsense
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
}*/

uint sliceSteps(float3 positionVS, float3 V, float2 start, float2 rayDir, float t, float step, float samplingDirection, float N, uint bitfield) {
    for (uint i = 0; i < STEPS; i++, t += step) {
        float2 samplePos = clamp(start + t * rayDir, 1, BUFFER_SCREEN_SIZE - 2);
        samplePos = round(samplePos) + 0.5;
        float3 samplePosVS = zfw::uvToView(samplePos.xy / BUFFER_SCREEN_SIZE);
        float3 delta = samplePosVS - positionVS;
		// to testers, if you replace the following line with the one right after, you need to also adjust thickness accordingly. 
		//float2 fb = acos(float2(dot(normalize(delta), V), dot(normalize(delta - V * THICKNESS * length(samplePosVS)), V)));
		
	    float2 fb = acos(float2(dot(normalize(delta), V), dot(normalize(delta - V * THICKNESS), V)));
	    fb = saturate(((samplingDirection * -fb) - N + PI/2) / PI);
	    fb = samplingDirection >= 0 ? fb.yx : fb.xy;

   	 uint a = round(fb.x * SECTORS);
    	uint b = round((fb.y - fb.x) * SECTORS);
    	bitfield |= ((1 << b) - 1) << a;
    }
    return bitfield;
}


float gtao(float2 uv, float2 vpos) {
	float2 random = stbn(vpos);
	
	float2 start = vpos * 2.0;
	float3 positionVS = zfw::uvToView(uv);
	positionVS.z *= 0.9999; // Move center pixel towards camera a bit.

	float ao = 0.0;
	
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
		aoBF = sliceSteps(positionVS, V, start, direction, offset, step, 1, N, aoBF);
		aoBF = sliceSteps(positionVS, V, start, -direction, offset, step, -1, N, aoBF);

		ao += float(countbits(aoBF));
	}
	ao = 1.0 - ao / (float(SECTORS) * SLICES);
	return positionVS.z > FAR_CLIP || ao < -0.001 ? 1.0 : ao;
}

float4 main(float4 vpos : SV_Position, float2 uv : TEXCOORD) : SV_Target {
	//float3 MV = zfw::getVelocity(uv);
	//float weight = 1.0 / tex2Dfetch(sAccum, vpos.xy); //bitch you will forget to div here and in read + save!!!!
	//return gtao(uv, vpos.xy) + tex2D(AOSPrev, uv + MV.xy).x * weight * MV.z; not yet
	//return lerp(gtao(uv, vpos.xy), tex2D(AOSPrev, uv + MV.xy).x, historySize * MV.z);
	float ao = gtao(uv, vpos.xy);
	return float4(ao, ao, ao, 1.0);
}

/*
uint read(float4 vpos : SV_Position, float2 uv : TEXCOORD) : SV_Target {
	return tex2Dfetch(sAccum, vpos.xy) + 1u;
}

uint save(float4 vpos : SV_Position, float2 uv : TEXCOORD) : SV_Target {
	return tex2Dfetch(sAccumW, vpos.xy);
}

float denoise0(float4 vpos : SV_Position, float2 uv : TEXCOORD) : SV_Target {
	return atrous(AOS, uv, 0);
}

float denoise1(float4 vpos : SV_Position, float2 uv : TEXCOORD) : SV_Target {
	return atrous(sDN0, uv, 1);
}

#define RCP_GAMMA 2.2
#define GAMMA 0.4545
float4 display(float4 vpos : SV_Position, float2 uv : TEXCOORD) : SV_Target {
	//float denoised = tex2Dfetch(sDN1, vpos.xy * 0.5).x;
	float denoised = JointBilateralUpsample(AOS, lowN, highN, uv).x;
	float4 bb = float4(zfw::getBackBuffer(uv), 1.0);
	bb = pow(bb, RCP_GAMMA);
	float luma = dot(bb, float4(0.2126, 0.7152, 0.0722, 0.));
	float protection = smoothstep(0.0, 1.0, luma);
	//denoised = lerp(denoised, 1.0, protection); // warning, your balls WILL explode!
	return pow(lerp(bb, lerp(denoised * bb, denoised, debug), strength), GAMMA);
}

float4 cache(float4 vpos : SV_Position, float2 uv : TEXCOORD) : SV_Target {
	return tex2Dfetch(AOS, vpos.xy).x;
}*/

technique SCVBAO {
	/*
	pass Read {
		VertexShader = PostProcessVS;
		PixelShader = read;
		RenderTarget = AccumW;
	}
	pass Increment {
		VertexShader = PostProcessVS;
		PixelShader = save;
		RenderTarget = Accum;
	}
	pass Denoise0 {
		VertexShader = PostProcessVS;
		PixelShader = denoise0;
		RenderTarget = DN0;
	}
	pass Denoise1 {
		VertexShader = PostProcessVS;
		PixelShader = denoise1;
		RenderTarget = DN1;
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
	}*/
	pass Main { 
		VertexShader = PostProcessVS;
		PixelShader = main;
	}
}
