#include "ReShade.fxh"
#include "soupcan_includes/ao.fxh"

float4 main(float4 vpos : SV_Position, float2 texcoord : TEXCOORD) : SV_Target {
	//return ReShade::GetLinearizedDepth(texcoord);
	return occlusion(texcoord, vpos.xy);// * tex2Dfetch(ReShade::BackBuffer, vpos.xy);
}

technique SoupCanAO {
	pass AO {
		VertexShader = PostProcessVS;
		PixelShader = main;
	}
}