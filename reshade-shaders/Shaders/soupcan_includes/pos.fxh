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

float3 GetScreenSpaceNormal(float2 texcoord) {
    float3 offset = float3(BUFFER_PIXEL_SIZE, 0.0);
    float2 posCenter = texcoord.xy;
    float2 posNorth  = posCenter - offset.zy;
    float2 posEast   = posCenter + offset.xz;

    float3 vertCenter = float3(posCenter - 0.5, 1) * ReShade::GetLinearizedDepth(posCenter);
    float3 vertNorth  = float3(posNorth - 0.5,  1) * ReShade::GetLinearizedDepth(posNorth);
    float3 vertEast   = float3(posEast - 0.5,   1) * ReShade::GetLinearizedDepth(posEast);

    return normalize(cross(vertCenter - vertNorth, vertCenter - vertEast));
}

float3 GetWorldSpaceNormal(float2 texcoord) {
    float3 offset = float3(BUFFER_PIXEL_SIZE, 0.0);
    float2 posCenter = texcoord.xy;
    float2 posNorth  = posCenter - offset.zy;
    float2 posEast   = posCenter + offset.xz;
	float3 vertCenter = getWorldPosition(posCenter,       GetDepth(posCenter));
    float3 vertNorth  = getWorldPosition(posNorth,  GetDepth(posNorth));
    float3 vertEast   = getWorldPosition(posEast,   GetDepth(posEast));
    
    return normalize(cross(vertCenter - vertNorth, vertCenter - vertEast));
}