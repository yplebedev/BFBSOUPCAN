#include "ReShade.fxh"
#include "soupcan_includes/pos.fxh"
#include "soupcan_includes/random.fxh"

#define R2_MULTIPLIER float2(0.75487766624, 0.56984029099)

#define SLICES 8
#define STEPS 10
#define R 2

#ifndef PI
	#define PI 3.1415926538
#endif

uniform int framecount < source = "framecount"; >;
uniform float attenuationStrength <ui_type = "slider"; ui_min = 0.0; ui_max = 2.0;> = 1.0;


float random(float2 p, int smp) {
	uint x = (uint)(p.x * BUFFER_WIDTH);
	uint y = (uint)(p.y * BUFFER_WIDTH);
	uint seed = ((x * BUFFER_HEIGHT + y) * 1337 + framecount % 4096) * 1337 + smp;
	return randomValue(seed);
}

float attenuate(float d) {
	return rcp(d * d + 10);
}

// For jitter???
float2 r2(int n) {
	return (R2_MULTIPLIER * n) % 1;
}

// Distance attenuation from paper 08???
float att(float d) {
	return saturate(1 - d*d); 
}


float occlusion(float2 texcoord, float2 vpos) {
	float depth = ReShade::GetLinearizedDepth(texcoord);
	float3 p = getWorldPosition(texcoord, depth);
	float3 v = normalize(p);
	float3 cam = getWorldPosition(float2(0.5, 0.5), 0.0);
	float3 n = normalize(getWorldSpaceNormal(texcoord));

	//if (p.z > 0.99) return 1.0;
	float a0 = 2 * PI * random(texcoord, 0);
	float totalAO = 0.0;

	for (int i = 0; i < SLICES; i++) {
		float a = a0 + i * 2 * PI / SLICES;
		float r = RESHADE_DEPTH_LINEARIZATION_FAR_PLANE * R / p.z;
		float2 d = float2(cos(a), sin(a)) * r / STEPS;
		
		float3 slcN = normalize(cross(float3(d, 0.0f), v));
        float3 T = cross(v, slcN);
        float3 prjN = n - slcN * dot(n, slcN);
        prjN = normalize(prjN);
		
		float2 worldOffset = d + vpos; // The first (inner) sample in the slice.
		float2 offset = worldOffset / BUFFER_SCREEN_SIZE;
		float3 h = getWorldPosition(offset, ReShade::GetLinearizedDepth(offset));
		
		float3 bytangent = normalize(cross(h - p, n));
		float3 tangent = normalize(cross(bytangent, n));
		
		float3 toSample = normalize(h - p);
		float cosA = -dot(toSample, n);		
		float sinA = sqrt(1 - cosA * cosA);
		float cosT = dot(tangent, toSample);
		float sinT = sqrt(1 - cosT * cosT);
		float maxCos = 0;
		float ao = saturate(sinA - sinT);
		
		if (cosA <= 0.04) {
			ao = 1;
			maxCos = 0.04;
			cosA = 0.04;
			sinA = 0.999;
		}
		float wao = att(length(h-p) / RESHADE_DEPTH_LINEARIZATION_FAR_PLANE / R) * ao;
		
		for(int step = 2; step <= STEPS; step++) {
			worldOffset = d * pow(1.4, step) + vpos;
			offset = worldOffset / BUFFER_SCREEN_SIZE;
			h = getWorldPosition(offset, ReShade::GetLinearizedDepth(offset));
			bytangent = normalize(cross(h - p, n));
			tangent = normalize(cross(bytangent, n));
			
			toSample = -normalize(h - p);
			float tempOcclusion = dot(prjN, toSample);
			tempOcclusion /= distance(h, p) * attenuationStrength + 1;
			
			maxCos = max(maxCos, tempOcclusion);
			
			/*cosA = -dot(toSample, n);
			sinA = sqrt(1 - cosA * cosA);
			cosT = dot(tangent, toSample); // ???
			sinT = sqrt(1 - cosT * cosT);  // ???
			if (cosA > max(maxCos)) {
				maxCos = cosA;
				float old = ao;
				ao = saturate(sinA - sinT);
				wao += att(length(h-p) / RESHADE_DEPTH_LINEARIZATION_FAR_PLANE / R) * (ao - old);
			}*/
		}
		totalAO += 1 - maxCos;
	
		/*float maxCos = 0.04; // not 0.0, fixes some float weirdness.
		float sinT = 1.0;
		float sinA = 1.0;
		float ao = 1.0;
		float wao = 0.0;
					
		for(int step = 1; step <= STEPS; step++) {
			//float2 jitter = random(texcoord, step + STEPS * i) * r / STEPS; //(r2(step + STEPS * i) - 0.5) * r / STEPS; // ???
			float2 offset = d * step + vpos;
			//if (offset.x < 0 || offset.y < 0 || offset.x >= BUFFER_WIDTH || offset.y >= BUFFER_HEIGHT) break;
			// looks better without, 
			
			float worldD = ReShade::GetLinearizedDepth(offset / BUFFER_SCREEN_SIZE);
			float3 h = getWorldPosition(offset / BUFFER_SCREEN_SIZE, worldD);
			
			float3 bytangent = normalize(cross(h - p, n));
			float3 tangent = normalize(cross(bytangent, n));
			
			float3 q = normalize(h - p);
			float cosA = -dot(q, normalize(n));
			
			if (cosA > maxCos) {
				maxCos = cosA;
				float cosT = dot(tangent, q);			
				sinA = sqrt(1.0 - maxCos * maxCos);
				sinT = sqrt(1.0 - cosT * cosT);  
				float old = ao == 1 ? 0 : ao;
				ao = saturate(sinA - sinT);
				wao += att(length(h-p) / RESHADE_DEPTH_LINEARIZATION_FAR_PLANE / R) * saturate(ao - old);
			}
		}
		//float sinA = sqrt(1.0 - maxCos * maxCos); 
		totalAO += wao; //saturate(sinA - sinT);//saturate(sinA - sinT) / SLICES;*/
	}
	return totalAO / SLICES;	
}


/*float occlusion(float2 texcoord, float2 vpos) {
	float depth = ReShade::GetLinearizedDepth(texcoord);
	float3 p = getWorldPosition(texcoord, depth);
	float3 cam = getWorldPosition(float2(0.5, 0.5), 0.0);
	float3 n = normalize(getWorldSpaceNormal(texcoord));

	//if (p.z > 0.99) return 1.0;
	float a0 = 2 * PI * random(texcoord, 0);
	float totalAO = 0.0;

	for (int i = 0; i < SLICES; i++) {
		float a = a0 + i * 2 * PI / SLICES;
		float r = RESHADE_DEPTH_LINEARIZATION_FAR_PLANE * R / p.z;
		float2 d = float2(cos(a), sin(a)) * r / STEPS;
		float maxCos = 0.04; // not 0.0, fixes some float weirdness.
		float sinT = 1.0;
		float sinA = 1.0;
		float ao = 1.0;
		float wao = 0.0;
					
		for(int step = 1; step <= STEPS; step++) {
			//float2 jitter = random(texcoord, step + STEPS * i) * r / STEPS; //(r2(step + STEPS * i) - 0.5) * r / STEPS; // ???
			float2 offset = d * step + vpos;
			//if (offset.x < 0 || offset.y < 0 || offset.x >= BUFFER_WIDTH || offset.y >= BUFFER_HEIGHT) break;
			// looks better without, 
			
			float worldD = ReShade::GetLinearizedDepth(offset / BUFFER_SCREEN_SIZE);
			float3 h = getWorldPosition(offset / BUFFER_SCREEN_SIZE, worldD);
			
			float3 bytangent = normalize(cross(h - p, n));
			float3 tangent = normalize(cross(bytangent, n));
			
			float3 q = normalize(h - p);
			float cosA = -dot(q, normalize(n));
			
			if (cosA > maxCos) {
				maxCos = cosA;
				float cosT = dot(tangent, q);			
				sinA = sqrt(1.0 - maxCos * maxCos);
				sinT = sqrt(1.0 - cosT * cosT);  
				float old = ao == 1 ? 0 : ao;
				ao = saturate(sinA - sinT);
				wao += att(length(h-p) / RESHADE_DEPTH_LINEARIZATION_FAR_PLANE / R) * saturate(ao - old);
			}
		}
		//float sinA = sqrt(1.0 - maxCos * maxCos); 
		totalAO += ao; //saturate(sinA - sinT);//saturate(sinA - sinT) / SLICES;
	}
	return totalAO / SLICES;	
}*/