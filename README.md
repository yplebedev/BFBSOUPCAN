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
  
SoupCanTraceRTGI:
  WARNING: NOT EVEN REMOTELY FINISHED. HERE FOR CODE HISTORY AND PLAYTESTING. Updates in the discord server!

![cachedImage](https://github.com/user-attachments/assets/81fd0891-934e-40ea-81c0-22027f4378e3)


