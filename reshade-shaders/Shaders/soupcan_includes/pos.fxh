#define BUFFER_SIZE3 int3(BUFFER_WIDTH,BUFFER_HEIGHT,RESHADE_DEPTH_LINEARIZATION_FAR_PLANE)

float3 getWorldPosition(float2 coords,float depth) {
        float3 result = float3((coords-0.5)*depth,depth);

        result *= BUFFER_SIZE3;
        return result;
}

float3 getScreenPosition(float3 world) {
	float3 result = world / BUFFER_SIZE3;
	return float3(result.xy/result.z + 0.5, result.z);
}