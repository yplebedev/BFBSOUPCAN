#define BUFFER_SIZE3 int3(BUFFER_WIDTH,BUFFER_HEIGHT,RESHADE_DEPTH_LINEARIZATION_FAR_PLANE)
#define GetDepth ReShade::GetLinearizedDepth

float3 getWorldPosition(float2 coords, float depth) {
        float3 result = float3((coords - 0.5) * depth, depth);

        result *= int3(BUFFER_WIDTH,BUFFER_HEIGHT,RESHADE_DEPTH_LINEARIZATION_FAR_PLANE);
        return result;
}

float3 getScreenPosition(float3 world) {
	float3 result = world / int3(BUFFER_WIDTH,BUFFER_HEIGHT,RESHADE_DEPTH_LINEARIZATION_FAR_PLANE);
	return float3(result.xy/result.z + 0.5, result.z);
}


float3 getScreenSpaceNormal(float2 texcoord) {
    float3 offset = float3(BUFFER_PIXEL_SIZE, 0.0);
    float2 posCenter = texcoord.xy;
    float2 posNorth  = posCenter - offset.zy;
    float2 posEast   = posCenter + offset.xz;

    float3 vertCenter = float3(posCenter - 0.5, 1) * ReShade::GetLinearizedDepth(posCenter);
    float3 vertNorth  = float3(posNorth - 0.5,  1) * ReShade::GetLinearizedDepth(posNorth);
    float3 vertEast   = float3(posEast - 0.5,   1) * ReShade::GetLinearizedDepth(posEast);

    return normalize(cross(vertNorth - vertCenter, vertCenter - vertEast));
}

float3 getWorldSpaceNormal(float2 texcoord) {
    float3 offset = float3(BUFFER_PIXEL_SIZE, 0);
    float2 posCenter = texcoord.xy;
    float2 posNorth  = posCenter - offset.zy;
    float2 posEast   = posCenter + offset.xz;
	float2 posSouth  = posCenter + offset.zy;
	float2 posWest   = posCenter - offset.xz;
	float  depthCenter = GetDepth(posCenter);
	float  depthNorth = GetDepth(posNorth);
	float  depthEast = GetDepth(posEast);
	float  depthSouth = GetDepth(posSouth);
	float  depthWest = GetDepth(posWest);
	bool2 dir = bool2(abs(depthCenter - depthNorth) < abs(depthCenter - depthSouth), abs(depthCenter - depthEast) < abs(depthCenter - depthWest));
	
	float3 vertCenter = getWorldPosition(posCenter,       depthCenter);
	/*
float3 vertNorth  = getWorldPosition(posNorth,  GetDepth(posNorth));
    	float3 vertEast   = getWorldPosition(posEast,   GetDepth(posEast));
    	if (true) return normalize(cross(vertCenter - vertNorth, vertCenter - vertEast));
*/
if (dir.x && dir.y) {
		float3 vertNorth  = getWorldPosition(posNorth,  GetDepth(posNorth));
    	float3 vertEast   = getWorldPosition(posEast,   GetDepth(posEast));
    	return -normalize(cross(vertCenter - vertNorth, vertCenter - vertEast));
	} else if (!dir.x && !dir.y) {
		float3 vertSouth  = getWorldPosition(posSouth,  GetDepth(posSouth));
    	float3 vertWest   = getWorldPosition(posWest,   GetDepth(posWest));
    	return -normalize(cross(vertCenter - vertSouth, vertCenter - vertWest));		
	} else if (dir.x) {
		float3 vertNorth  = getWorldPosition(posNorth,  GetDepth(posNorth));
		float3 vertWest   = getWorldPosition(posWest,   GetDepth(posWest));
		return -normalize(-cross(vertCenter - vertNorth, vertCenter - vertWest));
	} else {
		float3 vertSouth  = getWorldPosition(posSouth,  GetDepth(posSouth));
		float3 vertEast   = getWorldPosition(posEast,   GetDepth(posEast));
		return -normalize(-cross(vertCenter - vertSouth, vertCenter - vertEast));
	}
}

float4 getPackedWorldSpaceNormal(float2 texcoord) {
    float3 offset = float3(BUFFER_PIXEL_SIZE, 0);
    float2 posCenter = texcoord.xy;
    float2 posNorth  = posCenter - offset.zy;
    float2 posEast   = posCenter + offset.xz;
	float2 posSouth  = posCenter + offset.zy;
	float2 posWest   = posCenter - offset.xz;
	float  depthCenter = GetDepth(posCenter);
	float  depthNorth = GetDepth(posNorth);
	float  depthEast = GetDepth(posEast);
	float  depthSouth = GetDepth(posSouth);
	float  depthWest = GetDepth(posWest);
	bool2 dir = bool2(abs(depthCenter - depthNorth) < abs(depthCenter - depthSouth), abs(depthCenter - depthEast) < abs(depthCenter - depthWest));
	
	float3 vertCenter = getWorldPosition(posCenter,       depthCenter);
	/*
float3 vertNorth  = getWorldPosition(posNorth,  GetDepth(posNorth));
    	float3 vertEast   = getWorldPosition(posEast,   GetDepth(posEast));
    	if (true) return normalize(cross(vertCenter - vertNorth, vertCenter - vertEast));
*/
if (dir.x && dir.y) {
		float3 vertNorth  = getWorldPosition(posNorth,  GetDepth(posNorth));
    	float3 vertEast   = getWorldPosition(posEast,   GetDepth(posEast));
    	return float4(-normalize(cross(vertCenter - vertNorth, vertCenter - vertEast)), depthCenter);
	} else if (!dir.x && !dir.y) {
		float3 vertSouth  = getWorldPosition(posSouth,  GetDepth(posSouth));
    	float3 vertWest   = getWorldPosition(posWest,   GetDepth(posWest));
    	return float4(-normalize(cross(vertCenter - vertSouth, vertCenter - vertWest)), depthCenter);		
	} else if (dir.x) {
		float3 vertNorth  = getWorldPosition(posNorth,  GetDepth(posNorth));
		float3 vertWest   = getWorldPosition(posWest,   GetDepth(posWest));
		return float4(-normalize(-cross(vertCenter - vertNorth, vertCenter - vertWest)), depthCenter);
	} else {
		float3 vertSouth  = getWorldPosition(posSouth,  GetDepth(posSouth));
		float3 vertEast   = getWorldPosition(posEast,   GetDepth(posEast));
		return float4(-normalize(-cross(vertCenter - vertSouth, vertCenter - vertEast)), depthCenter);
	}
}

float3 getWorldSpaceNormal2(float2 texcoord) {
    float3 offset = float3(BUFFER_PIXEL_SIZE, 0.0);
    float2 posCenter = texcoord.xy;
    float2 posNorth  = posCenter + offset.zy;
    float2 posEast   = posCenter - offset.xz;
	float3 vertCenter = getWorldPosition(posCenter,       GetDepth(posCenter));
    float3 vertNorth  = getWorldPosition(posNorth,  GetDepth(posNorth));
    float3 vertEast   = getWorldPosition(posEast,   GetDepth(posEast));
    
    return normalize(cross(vertCenter - vertNorth, vertCenter - vertEast));
}