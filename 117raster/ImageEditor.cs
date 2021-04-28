// Lukáš Salak
// Default = color offset correction function.

formula.contextCreate = (in Bitmap input, in string param) => {
if (string.IsNullOrEmpty(param))
  return null;

Dictionary<string, string> p = Util.ParseKeyValueList(param);
Dictionary<string, object> sc = new Dictionary<string, object>();


float Brightness = 0.0f;
if (Util.TryParse(p, "Brightness", ref Brightness))
  Brightness = Util.Clamp(Brightness, -100.0f, 100.0f);

float Contrast = 0.0f;
if (Util.TryParse(p, "Contrast", ref Contrast))
  Contrast = Util.Clamp(Contrast, -100.0f, 100.0f);

float Highlights = 0.0f;
if (Util.TryParse(p, "Highlights", ref Highlights))
  Highlights = Util.Clamp(Highlights, -100.0f, 100.0f);

float Shadows = 0.0f;
if (Util.TryParse(p, "Shadows", ref Shadows))
  Shadows = Util.Clamp(Shadows, -100.0f, 100.0f);

float Saturation = 0.0f;
if (Util.TryParse(p, "Saturation", ref Saturation))
  Saturation = Util.Clamp(Saturation, -100.0f, 100.0f);

sc["Brightness"] = Brightness;
sc["Contrast"] = Contrast;
sc["Highlights"] = Highlights;
sc["Shadows"] = Shadows;
sc["Saturation"] = Saturation;
sc["tooltip"] = "Brightness=<int> |  Changes the brightness of the image (-100 to 100)\r"+
				"Contrast=<int> |  Changes the contrast of the image (-100 to 100)\r"+
				"Highlights=<int> |  Changes the white balance of the image (-100 to 100)\r"+
				"Shadows=<int> |  Changes the black balance of the image (-100 to 100)\r"+
				"Saturation=<int> |  Changes the saturation of the image (-100 to 100)\r";


return sc;
};


formula.pixelTransform0 = (
in ImageContext ic,
ref float R,
ref float G,
ref float B) =>
{
	
List<float> channels = new List<float>(){R,G,B};

float luminance = (float)Math.Sqrt( Math.Pow(R,2.0f) + Math.Pow(G,2.0f) + Math.Pow(B,2.0f));


for(int i = 0; i < 3; i++){
	
// Brightness	
float BrightnessRatio = (float)ic.context["Brightness"]/100.0f;
// Contrast
float ContrastRatio = ((259.0f/255.0f)*(( (float)ic.context["Contrast"]/100.0f) +1)/((259.0f/255.0f)-(float)ic.context["Contrast"]/100.0f) ) * (channels[i] - 0.5f) + 0.5f;	
// Highlights
float high = ((float)ic.context["Highlights"]/100.0f) * 0.05f * ( (float)Math.Pow(8.0f, luminance) - 1.0f );
// Shadows
float shad = ((float)ic.context["Shadows"]/100.0f) * 0.05f * ( (float)Math.Pow(8.0f, 1.0f - luminance) - 1.0f );

	channels[i] = ContrastRatio + BrightnessRatio + high - shad;	
}



R = channels[0];
G = channels[1]; 
B = channels[2];

//Saturation
Arith.RGBtoHSV(R,G,B, out float Hue, out float Saturation, out float Value);
Saturation *= (float)Math.Pow(1 + (float)ic.context["Saturation"]/100.0f, 2);
Arith.HSVtoRGB(Hue,Saturation,Value, out R, out G, out B);


return true;
};
