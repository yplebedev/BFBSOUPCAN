#include "ReShade.fxh"
#include "soupcan_includes/pos.fxh"
#include "soupcan_includes/random.fxh"

#define getDepth ReShade::GetLinearizedDepth

#define SLICES 8
#define STEPS 5
#define R 4
#define R_MAX_CLAMP 1200
#define THICKNESS 0.3
#define SECTOR_COUNT 32 // chaning this makes the error "switch". try 16 and 32.
#define FAR_CLIP (RESHADE_DEPTH_LINEARIZATION_FAR_PLANE - 1)


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

// Modified from https://www.shadertoy.com/view/4djSRW
// It's called once per pixel.
float2 fastHash(float2 p) {
	float3 p3 = float3(p * BUFFER_SCREEN_SIZE, framecount & 0xFFFF);
	p3 = frac(p3 * float3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx + 33.33);
    return frac((p3.xy + p3.yz) * p3.zx);
}

uint sliceSteps(float3 positionVS, float3 V, float2 start, float2 rayDir, float t, float step, float samplingDirection, float N, uint bitfield) {
	for (uint i = 0; i < STEPS; i++, t += step) {
		float2 samplePos = clamp(start + t * rayDir, 1, BUFFER_SCREEN_SIZE - 2);
		float3 samplePosVS = getWorldPosition(samplePos.xy / BUFFER_SCREEN_SIZE, getDepth(samplePos));
		float3 delta = samplePosVS - positionVS;
		
		float2 fb = acos(float2(dot(normalize(delta), V), dot(normalize(delta - V * THICKNESS), V)));
		fb = saturate(((samplingDirection * -fb) - N + PI/2) / PI);
		fb = samplingDirection >= 0 ? fb.yx : fb.xy;
	
		uint a = fb.x * SECTOR_COUNT;
		uint b = ceil((fb.y - fb.x) * SECTOR_COUNT);
		bitfield |= b > 0 ? (0xFFFFFFFF >> (SECTOR_COUNT - b)) << a : 0;
	}
		return bitfield;
}

//float stepBaseNoise = (slice + step * STEPS) * 0.6180339887498948482; 	
//float stepNoise = 2 * frac(tex2Dfetch(bn_sampler, vpos % 64).y + stepBaseNoise);

float vbgtao(float2 uv, float2 vpos) {
	float2 random = fastHash(uv);

	float depth = GetDepth(uv);
	float3 cPosV = getWorldPosition(uv, depth);
	cPosV.z *= 0.99999;  // Move center pixel slightly towards camera
	// to avoid imprecision artifacts due to depth buffer imprecision.
	float3 viewV = normalize(-cPosV);
	float3 normalV = getWorldSpaceNormal(uv);
	float visibility = 0;
	
	float randomAngle = tex2Dfetch(bn_sampler, vpos % 64).x;
	
	float scaling = R / length(cPosV);
	scaling = clamp(scaling, 0, R);
	
	for (float slice = 0; slice < SLICES; slice++) {
		float phi = PI * (slice / SLICES + randomAngle);

		float2 omega = float2(cos(phi), sin(phi));
		//omega.y = omega.y * BUFFER_ASPECT_RATIO;
		
		float3 dirV = float3(omega, 0);
		float3 orthoDirV = dirV - dot(dirV, viewV) * viewV;
		float3 axisV = cross(dirV, viewV);
		float3 projN = normalV - axisV * dot(normalV, axisV);
		
		int sgnN = sign(dot(orthoDirV, projN));
		float cosN = saturate(dot(projN, viewV) / length(projN));
		float n = sgnN * acos(cosN);
		
		uint bitmask = 0u;
		float offset = max(random.y * scaling, length(BUFFER_PIXEL_SIZE));
		bitmask = sliceSteps(cPosV, viewV, vpos, omega, offset, scaling, 1, n, bitmask);
		bitmask = sliceSteps(cPosV, viewV, vpos, -omega, offset, scaling, -1, n, bitmask);
		
		visibility += float(countbits(bitmask)) / float(SECTOR_COUNT);
	}
	visibility /= SLICES;
	return cPosV.z > FAR_CLIP || visibility < -0.001 ? 1.0 : visibility;
}

float vbgtao_other(float2 uv) {
	float2 random = fastHash(uv);

	float2 start = uv * BUFFER_SCREEN_SIZE; // this is vpos
	float3 positionVS = getWorldPosition(uv);
	positionVS.z *= 0.9999; // Move center pixel towards camera a bit.

	float ao = 0.0;
	
	float3 V = normalize(-positionVS);
	float3 normalVS = getWorldSpaceNormal(uv);

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

		ao += 1.0 - float(countbits(aoBF)) / float(SECTOR_COUNT);
	}
	ao /= SLICES;
	return positionVS.z > FAR_CLIP || ao < -0.001 ? 1.0 : ao;
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
				
				currentCos = lerp(currentCos, 0, saturate(distance(sPosV, cPosV) / 10));
				
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