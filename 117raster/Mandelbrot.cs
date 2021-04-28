// Text params -> script context.
// Any global pre-processing is allowed here.

// Lukáš Salak

formula.contextCreate = (in Bitmap input, in string param) => {
  if (string.IsNullOrEmpty(param))
    return null;

  Dictionary<string, string> p = Util.ParseKeyValueList(param);
  Dictionary<string, object> sc = new Dictionary<string, object>();


  bool fill = false;
  Util.TryParse(p, "fill", ref fill);

  bool rainbow = false;
  Util.TryParse(p, "rainbow", ref rainbow);


  int depth = 10;
  if (Util.TryParse(p, "depth", ref depth))
    depth = Util.Clamp(depth, 1, 1000);

  float limit = 10;
  if (Util.TryParse(p, "limit", ref limit))
    limit = Util.Clamp(limit, 1, 1000);


  float xBoundary = 2.0f;
  if (Util.TryParse(p, "xBoundary", ref xBoundary))
    xBoundary = Util.Clamp(xBoundary, 0.1f, 10.0f);

  float yBoundary = 2.0f;
  if (Util.TryParse(p, "yBoundary", ref yBoundary))
    yBoundary = Util.Clamp(yBoundary, 0.1f, 10.0f);

  sc["limit"] = limit;
  sc["rainbow"] = rainbow;
  sc["depth"] = depth;
  sc["fill"] = fill;
  sc["xBoundary"] = xBoundary;
  sc["yBoundary"] = yBoundary;
  sc["tooltip"] = "rainbow=<bool> .. sets rainbow effect (true or false)\r" +
                  "fill=<bool> .. fill the object or not (what dont go over the limit) (true or false)\r" +
                  "limit=<float> .. limit in function (what is more than limit) (1 to 200)\r" +
                  "depth=<int> .. depth of calculation (1 to 1000)\r" +
                  "xBoundary=<float> .. range of calculation in x axis (0.1 to 10)\r" +
                  "yBoundary=<float> .. range of calculation in y axis (0.1 to 10)\r";

  return sc;
};


// Test create function: mandelbrot set
formula.pixelCreate = (
  in ImageContext ic,
  out float R,
  out float G,
  out float B) => {

    ic.context.TryGetValue("depth", out object depth);
    ic.context.TryGetValue("limit", out object limit);

    ic.context.TryGetValue("fill", out object fill);
    ic.context.TryGetValue("rainbow", out object rainbow);

    ic.context.TryGetValue("xBoundary", out object xBoundary);
    ic.context.TryGetValue("yBoundary", out object yBoundary);
    float normalizedX = -(float)xBoundary + ic.x * (((float)xBoundary + (float)xBoundary - 1.2f / (float)xBoundary) / ic.width);
    float normalizedY = -(float)yBoundary + ic.y * (((float)yBoundary + (float)yBoundary) / ic.height);

    var baseX = normalizedX;
    var baseY = normalizedY;
    float realsquared;
    float imaginarysquared;

    int i = 0;
    for (; i < (int)depth; i++) {

      // (a+bi)*(a+bi) = a^2 - b^2 + 2ab
      realsquared = normalizedX * normalizedX - normalizedY * normalizedY;
      imaginarysquared = 2 * normalizedX * normalizedY;

      normalizedX = realsquared + baseX;
      normalizedY = imaginarysquared + baseY;

      if (System.Math.Abs(normalizedX + normalizedY) > (float)limit) break; // limit is my infinity here
    }


    float tempI = i;
    float tempDepth = (int)depth;
    float brightness = tempI / tempDepth;
    if (i == (int)depth && !(bool)fill) {
      brightness = 0.0f; }


    // could have done some more thought into doing this but i like this
    if ((bool)rainbow) {
      R = brightness * ((float)ic.x / (float)ic.width);
      G = brightness * ((float)ic.y / (float)ic.height);
      B = brightness * ((float)System.Math.Abs(ic.x - ic.y) / (float)(ic.width+ic.height));
    }
      else {
      R = brightness;
      G = brightness;
      B = brightness;
    }
};
