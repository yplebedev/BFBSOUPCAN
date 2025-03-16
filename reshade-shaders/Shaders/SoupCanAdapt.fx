#include "ReShade.fxh"
#include "soupcan_includes/blur.fxh"
#include "soupcan_includes/color_ops.fxh"

/*
#ifndef UjelHDR_PARITY
	#define UjelHDR_PARITY 0
#endif
some day...
*/ 

uniform int channel <ui_type = "combo";
					 ui_label = "Adaptation source";
					 ui_tooltip = "Where to source the brightness from. Choose what looks best.";
					 ui_items = "Red\0Green\0Blue\0Luminance\0";
					 ui_type = "slider"; > = 0;
					

uniform float strength <  ui_label = "Adaptation strength";
					  	ui_tooltip = "How strong the effect should be. Choose whatever doesn't flatten the image too much.";
					  	ui_type = "slider"; > = 0.2;
					

float4 multiplyMix(float4 under, float4 over, float alpha) {
	return lerp(under, under * over, alpha);
}

float getAdapt(float4 blurred) {
	switch (channel) {
		case (0): return blurred.r;
		case (1): return blurred.g;
		case (2): return blurred.b;
		case (3): return dot(blurred.rgb, lumaCoeff);
	}
	return blurred.g;
}

float4 main(float4 vpos : SV_Position, float2 texcoord : TEXCOORD) : SV_Target {
	float4 input = sRGBtoLinear(tex2Dfetch(ReShade::BackBuffer, vpos.xy));
	float4 OKLabIn = lrgb2oklab(input);
	
	float4 blurred = tex2Dfetch(BSam5, vpos.xy);
	float adaptChannel = getAdapt(blurred);
	float4 extraBright = pow(input, 0.5);
	
	float4 adapted = multiplyMix(extraBright, input, adaptChannel);
	
	float4 OKLabOut = lrgb2oklab(adapted);
	float adaptedLuminance = lerp(OKLabOut.x, OKLabIn.x, 1 - strength);
	adapted = float4(adaptedLuminance, OKLabIn.yz, input.w);
	
	return linearTosRGB(oklab2lrgb(adapted));
}

float4 prepare(float4 vpos : SV_Position, float2 texcoord : TEXCOORD) : SV_Target {
	float4 backBufferRes = tex2Dfetch(ReShade::BackBuffer, vpos.xy);
	backBufferRes = sRGBtoLinear(backBufferRes);
	
	return backBufferRes;
}

technique SoupCanAdapt {
	pass prepass {
		VertexShader = PostProcessVS;
		PixelShader = prepare;
		RenderTarget = inTex;
	}
	pass {
		VertexShader = PostProcessVS;
		PixelShader = DownSample0;
		RenderTarget = BTex4;
	}
	pass {
		VertexShader = PostProcessVS;
		PixelShader = DownSample1;
		RenderTarget = BTex3;
	}
	pass {
		VertexShader = PostProcessVS;
		PixelShader = DownSample2;
		RenderTarget = BTex2;
	}
	pass {
		VertexShader = PostProcessVS;
		PixelShader = DownSample3;
		RenderTarget = BTex1;
	}
	pass {
		VertexShader = PostProcessVS;
		PixelShader = DownSample4;
		RenderTarget = BTex0;
	}
	pass {
		VertexShader = PostProcessVS;
		PixelShader = UpSample0;
		RenderTarget = BTex1;
	}
	pass {
		VertexShader = PostProcessVS;
		PixelShader = UpSample1;
		RenderTarget = BTex2;
	}
	pass {
		VertexShader = PostProcessVS;
		PixelShader = UpSample2;
		RenderTarget = BTex3;
	}
	pass {
		VertexShader = PostProcessVS;
		PixelShader = UpSample3;
		RenderTarget = BTex4;
	}
	pass {
		VertexShader = PostProcessVS;
		PixelShader = UpSample4;
		RenderTarget = BTex5;
	}
	pass main {
		VertexShader = PostProcessVS;
		PixelShader = main;
	}
}