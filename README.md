More of my contributions to the ReShade project.
Not affiliated with smallbbsoop's stuff, but I am a fan of her work.

COMMUNITY DEACHRONYMS:
  None. YET.

SHADERS USED: Before/After, UjelHDR, SoupCanAdapt.
![207243~1](https://github.com/user-attachments/assets/4b7e4d0c-ed71-4361-8710-e2261a91a5ca)


SoupCanAdapt:
An adaptation shader that makes your dark-er screenshots less dark. TO BE CLEAR, this IS (debatably, gamma fusion would be more correct) exposure fusion, but NOT done in the same way as Marty's or Zenteon's is.
When using ReShade, most effects, one way or another, increase global contrast. Not only is this a problem with "tonemapping" and "HDR" shaders like mine, but with even less realized effects, like AO and SSGI. Naturally, his becomes a problem if not fixed. Lowering the contrast naïvely does fix the issue, but it lowers LOCAL contrast, and that is usually desired, if not outright intentionally added.
This serves as a potential fix, inspired by a traditional Photoshop method, used widely in one way or another.

HOW IT WORKS:

  Step 1: Convert to linear to make sure everything is mixed correctly
  
  Step 2: Very scary looking blur
  
  Step 3: Use a multiply mix to combine a muted version of the image with the original, based on a single blurred channel, or luminance.



HOW TO USE:

  This is a suggestion for your workflow.
  
  Step 1: Enable the shader, and put it after your DoFs and lightning / depth-dependant shaders.
  
  Step 2: Switch the mix source until it looks fitting. Green works fine, luma works whenever others might now, although this might depend on the scene
  
  Step 3: Tone it up or down, based on taste and desired result
  
  Step 4 (in here for newbies): Do not touch the blur, unless it is too small, and you see halos.

Shaders used: Framework, DH_UBER_MASK, ReVeil, VBAO, iMMERSE SMAA, FSR1_2X. Credit to U.K.N.
![SCVBAO](https://github.com/user-attachments/assets/e66ae8bc-3aa5-482b-81f2-5ba77432b9e3)
SCAO:

An ambient occlusion that objectivly doesn't quite rival alternatives in "objective" measures (read: is kinda mid), but it exposes multiple preproc settings for nailing certain looks, in the limited testing I've recieved, it's much better at detail preservation. It fits a reasonable performance target, and scales from mid-range to high-end hardware reasonably. This effect has been developed free from ad-hoc falloff parameters, and instead depends only and exclusivly on visibility bitmasking to handle thin geometry accuratly. 
HOW IT WORKS:

This will be quite brief as the implementation has a lot of details that will perhaps one day I will cover it on more detail.
This follows the standard slice formulation, splooshing the integration dimention, that speeds up integration, lowers the amount of wasted memory fetches and improves the manifestation or blue noise in the final result. The difference between this technique and GTAO primarily boils down to the fact, that GTAO only really bothers with two angles to represent and entire plane's worth of geometry. This (in GTAO) leads to overoccludion behind close-by objects. Visibility bitmasking instead stores a bitmask, a field of 32 directions, either occluded or not. This allows to assume thickness, and have each step occlude a handful of directions at once. The resulting occlusion is calculated via a countbits() call, and, while being much harder to work with and a notch slower, the end result behaves much more reasonably, in my opinion.
The final image is computed by denoising and upscaling, as the entire image is generated in halfres. A joint bilateral filter is used for both. 

HOW TO USE:
  Step 1: Enable the shader and ZenteonFramework. Put framework before SCAO.
  Step 2: Adjust radius.
  Step 3: Profit.



Image albums, in case I make more:
https://imgsli.com/MzkxNzAy
