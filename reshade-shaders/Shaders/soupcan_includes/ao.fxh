#include "ReShade.fxh"
#include "soupcan_includes/pos.fxh"
#include "soupcan_includes/random.fxh"

#define SLICES 8
#define STEPS 5
#define R 3.0
#define THICKNESS 0.3
#define SECTOR_COUNT 32
#define NOISE_TEX_NAME "QuarkBN_SOUPGI_TINY.png"

texture bn_tex < source = NOISE_TEX_NAME; > { Width = 64; Height = 64; Format = RGBA8; };
sampler bn_sampler { Texture = bn_tex; MagFilter = POINT; MinFilter = POINT;
	AddressU = REPEAT;
	AddressV = REPEAT;
	AddressW = REPEAT; };
	

texture normalTex { Width = BUFFER_WIDTH; Height = BUFFER_HEIGHT; Format = RGBA16F; };
sampler normalTexSampler { Texture = normalTex; };


#ifndef PI
	#define PI 3.1415926538
#endif

uniform int framecount < source = "framecount"; >;

float IGN(int2 xy){
	float3 conVr = float3(0.06711056, 0.00583715, 52.9829189);
	return frac(conVr.z * frac(dot(xy,conVr.xy)));
}

float white(float2 p, int smp) {
	uint x = (uint)(p.x * BUFFER_WIDTH);
	uint y = (uint)(p.y * BUFFER_WIDTH);
	uint seed = ((x * BUFFER_HEIGHT + y) * 1337 + framecount % 4096) * 1337 + smp;
	return randomValue(seed);
}

float attenuate(float d, float m) {
	float t = rcp(d * d * m + 1);
	return t;
}


float map(float value, float min1, float max1, float min2, float max2) {
	float perc = (value - min1) / (max1 - min1);
	return perc * (max2 - min2) + min2;
}


float3 vbgtao(float2 uv, float2 vpos) {
	float3 cPosV = getWorldPosition(uv, GetDepth(uv));
	cPosV.z *= 0.99999;  // Move center pixel slightly towards camera
	// to avoid imprecision artifacts due to depth buffer imprecision.
	float3 viewV = normalize(-cPosV);
	float3 normalV = getWorldSpaceNormal(uv);
	float visibility = 0;
	
	float randomAngle = tex2Dfetch(bn_sampler, vpos % 64).x;
	for (float slice = 0; slice < SLICES; slice++) {
		float phi = PI * (slice / SLICES + randomAngle);

		float2 omega = float2(cos(phi), sin(phi));
		omega.y = omega.y * BUFFER_ASPECT_RATIO;
		
		float3 dirV = float3(omega, 0);
		float3 orthoDirV = dirV - dot(dirV, viewV) * viewV;
		float3 axisV = cross(dirV, viewV);
		float3 projN = normalV - axisV * dot(normalV, axisV);
		
		int sgnN = sign(dot(orthoDirV, projN));
		float cosN = saturate(dot(projN, viewV) / length(projN));
		float n = sgnN * acos(cosN);
		
		[unroll]
		for (int side = 0; side < 2; side++) {
			float cHorizonCos = -1;
			for (float step = 0; step < STEPS; step++) {
				float stepBaseNoise = (slice + step * STEPS) * 0.6180339887498948482; 	
				float stepNoise = 2 * frac(tex2Dfetch(bn_sampler, vpos % 64).y + stepBaseNoise);
				float s = (step + stepNoise) / STEPS;
				s = pow(s, 2.0) + length(BUFFER_PIXEL_SIZE);
				float scaling = R / length(cPosV);
				scaling = clamp(scaling, 0, R);
				float2 sTexCoord = uv + (-1.0 + 2.0 * side) * s * scaling * omega;
				float3 sPosV = getWorldPosition(sTexCoord, GetDepth(sTexCoord));
				float3 sDelta = sPosV - cPosV;		
							
				
				
				float3 sHorizonV = normalize(sDelta);
				float currentCos = dot(sHorizonV, viewV);
				cHorizonCos = max(cHorizonCos, currentCos);
			}
			
			float h = n + clamp((-1.0 + 2.0 * side) * acos(cHorizonCos) - n, -PI/2, PI/2);
			visibility += length(projN) * (cosN + 2.0 * h * sin(n) - cos(2.0 * h - n)) / 4.0;
		}
	}
	return visibility / SLICES;
}


float gtao(float2 uv, float2 vpos) {
	float depth = GetDepth(uv);
	float3 cPosV = getWorldPosition(uv, depth);
	cPosV.z *= 0.99999;  // Move center pixel slightly towards camera
	// to avoid imprecision artifacts due to depth buffer imprecision.
	float3 viewV = normalize(-cPosV);
	float3 normalV = getWorldSpaceNormal(uv);
	float visibility = 0;
	
	float randomAngle = tex2Dfetch(bn_sampler, vpos % 64).x;
	for (float slice = 0; slice < SLICES; slice++) {
		float phi = PI * (slice / SLICES + randomAngle);

		float2 omega = float2(cos(phi), sin(phi));
		omega.y = omega.y * BUFFER_ASPECT_RATIO;
		
		float3 dirV = float3(omega, 0);
		float3 orthoDirV = dirV - dot(dirV, viewV) * viewV;
		float3 axisV = cross(dirV, viewV);
		float3 projN = normalV - axisV * dot(normalV, axisV);
		
		int sgnN = sign(dot(orthoDirV, projN));
		float cosN = saturate(dot(projN, viewV) / length(projN));
		float n = sgnN * acos(cosN);
		
		[unroll]
		for (int side = 0; side < 2; side++) {
			float cHorizonCos = -1;
			for (float step = 0; step < STEPS; step++) {
				float stepBaseNoise = (slice + step * STEPS) * 0.6180339887498948482; 	
				float stepNoise = 2 * frac(tex2Dfetch(bn_sampler, vpos % 64).y + stepBaseNoise);
				float s = (step + stepNoise) / STEPS;
				s = pow(s, 2.0) + length(BUFFER_PIXEL_SIZE);
				float scaling = R / length(cPosV);
				scaling = clamp(scaling, 0, R);
				float2 sTexCoord = uv + (-1.0 + 2.0 * side) * s * scaling * omega;
				float3 sPosV = getWorldPosition(sTexCoord, GetDepth(sTexCoord));
				float3 sDelta = sPosV - cPosV;
				float3 sHorizonV = normalize(sDelta);
				float currentCos = dot(sHorizonV, viewV);
				
				currentCos = lerp(currentCos, 0, saturate(length(sDelta) / 10));
				
				cHorizonCos = max(cHorizonCos, currentCos);
			}
			
			float h = n + clamp((-1.0 + 2.0 * side) * acos(cHorizonCos) - n, -PI/2, PI/2);
			visibility += length(projN) * (cosN + 2.0 * h * sin(n) - cos(2.0 * h - n)) / 4.0;
		}
	}
	return visibility / SLICES;
}

float4 hbil(float2 uv, float2 vpos) {
	float3 cPosV = getWorldPosition(uv, GetDepth(uv));
	cPosV.z *= 0.99999;  // Move center pixel slightly towards camera
	// to avoid imprecision artifacts due to depth buffer imprecision.
	float3 viewV = normalize(-cPosV);
	//float3 normalV = getWorldSpaceNormal(uv);
	float3 normalV = tex2D(normalTexSampler, uv).xyz;
	float visibility = 0;
	float3 light = 0;
	
	float randomAngle = tex2Dfetch(bn_sampler, vpos % 64).x;
	for (float slice = 0; slice < SLICES; slice++) {
		float phi = PI * (slice / SLICES + randomAngle);

		float2 omega = float2(cos(phi), sin(phi));
		omega.y = omega.y * BUFFER_ASPECT_RATIO;
		
		float3 dirV = float3(omega, 0);
		float3 orthoDirV = dirV - dot(dirV, viewV) * viewV;
		float3 axisV = cross(dirV, viewV);
		float3 projN = normalV - axisV * dot(normalV, axisV);
		
		int sgnN = sign(dot(orthoDirV, projN));
		float cosN = saturate(dot(projN, viewV) / length(projN));
		float n = sgnN * acos(cosN);
		
		[unroll]
		for (int side = 0; side < 2; side++) {
			float cHorizonCos = -1;
			for (float step = 0; step < STEPS; step++) {
				float stepBaseNoise = (slice + step * STEPS) * 0.6180339887498948482; 	
				float stepNoise = 2 * frac(tex2Dfetch(bn_sampler, vpos % 64).y + stepBaseNoise);
				float s = (step + stepNoise) / STEPS;
				s = pow(s, 2.0) + length(BUFFER_PIXEL_SIZE);
				float scaling = R / length(cPosV);
				scaling = clamp(scaling, 0, R);
				float2 sTexCoord = uv + (-1.0 + 2.0 * side) * s * scaling * omega;
				float3 sPosV = getWorldPosition(sTexCoord, GetDepth(sTexCoord));
				float3 sDelta = sPosV - cPosV;
				float3 sHorizonV = normalize(sDelta);
				float currentCos = dot(sHorizonV, viewV);
				cHorizonCos = max(cHorizonCos, currentCos);
				
				float3 col = tex2Dfetch(ReShade::BackBuffer, sTexCoord * BUFFER_SCREEN_SIZE).rgb;
				float contribution = max(cHorizonCos - currentCos, 0);
				float3 normal = tex2D(normalTexSampler, sTexCoord).rgb;
				light += col * contribution * saturate(dot(normal, -normalV)) * saturate(dot(viewV, normalV));
			}
			
			float h = n + clamp((-1.0 + 2.0 * side) * acos(cHorizonCos) - n, -PI/2, PI/2);
			visibility += length(projN) * (cosN + 2.0 * h * sin(n) - cos(2.0 * h - n)) / 4.0;
		}
	}
	return float4(light / SLICES, visibility / SLICES);
} 