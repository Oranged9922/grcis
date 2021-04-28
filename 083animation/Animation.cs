//using CircleCanvas;
//using MathSupport;
//using System;
//using System.Collections.Generic;
//using System.Diagnostics;
//using System.Drawing;
//using Utilities;

//namespace _083animation
//{
//  public class Animation
//  {
//    public static Ball b { get; set; }
//    public static Dictionary<string, string> p;
//    public static double[] yVals { get; set; }

//    /// <summary>
//    /// Form data initialization.
//    /// </summary>
//    /// <param name="name">Your first-name and last-name.</param>
//    /// <param name="wid">Initial image width in pixels.</param>
//    /// <param name="hei">Initial image height in pixels.</param>
//    /// <param name="from">Start time (t0)</param>
//    /// <param name="to">End time (for animation length normalization).</param>
//    /// <param name="fps">Frames-per-second.</param>
//    /// <param name="param">Optional text to initialize the form's text-field.</param>
//    /// <param name="tooltip">Optional tooltip = param help.</param>
//    public static void InitParams (out string name, out int wid, out int hei, out double from, out double to, out double fps, out string param, out string tooltip)
//    {
//      // {{

//      // Put your name here.
//      name = "Lukáš Salak";

//      // Frame size in pixels.
//      wid = 640;
//      hei = 640;

//      // Animation.
//      from = 0.0;
//      to = 10.0;
//      fps = 60.0;

//      // Form params.
//      param = "HLines=10,VLines=10,speed=5,seed=224,shadowSize=1.3,Xrotation=-0.5,Yrotation=-0.2,Zrotation=-0.2,colorful=1";
//      tooltip = "HLines=<int>, " +
//                "VLines=<int>, " +
//                "speed=<double>, " +
//                "seed=<long>, " +
//                "X/Y/Zrotation<float> (-1 to 1), " +
//                "colorful=<int> (0 or 1)";

//      // }}
//    }

//    /// <summary>
//    /// Global initialization. Called before each animation batch
//    /// or single-frame computation.
//    /// </summary>
//    /// <param name="width">Width of the future canvas in pixels.</param>
//    /// <param name="height">Height of the future canvas in pixels.</param>
//    /// <param name="start">Start time (t0)</param>
//    /// <param name="end">End time (for animation length normalization).</param>
//    /// <param name="fps">Required fps.</param>
//    /// <param name="param">Text parameter field from the form.</param>
//    public static void InitAnimation (int width, int height, double start, double end, double fps, string param)
//    {
//      // {{ TODO: put your init code here

//      // Parse parameters.

//      int HLines = 15;
//      int VLines = 15;
//      long seed = 4;
//      float shadowSize = 0.5f;

//      p = Util.ParseKeyValueList(param);
//      if (Animation.p.Count > 0)
//      {
//        // objects=<int>
//        if (Util.TryParse(Animation.p, "HLines", ref HLines))
//        {
//          if (HLines < 5)
//            HLines = 5;
//        }
//        if (Util.TryParse(Animation.p, "VLines", ref VLines))
//        {
//          if (VLines < 5)
//            VLines = 5;
//        }

//        if (Util.TryParse(Animation.p, "shadowSize", ref shadowSize))
//        {
//          if (shadowSize < 0.1f)
//            shadowSize = 0.1f;
//        }
//        Util.TryParse(Animation.p, "seed", ref seed);

//      }

//      Animation.b = new Ball(1, HLines, VLines, seed, shadowSize);
//      yVals = new double[b.Yvalues.Count];
//    }
//    // }}

//    /// <summary>
//    /// Draw single animation frame.
//    /// </summary>
//    /// <param name="c">Canvas to draw to.</param>
//    /// <param name="time">Current time in seconds.</param>
//    /// <param name="start">Start time (t0)</param>
//    /// <param name="end">End time (for animation length normalization).</param>
//    /// <param name="param">Optional string parameter from the form.</param>
//    public static void DrawFrame (Canvas c, double time, double start, double end, string param)
//    {
//      // {{ TODO: put your drawing code here

//      c.Clear(Color.Black);
//      c.SetAntiAlias(true);
//      float speed = 1;
//      float Xrotation = 0;
//      float Yrotation = 0;
//      float Zrotation = 0;

//      int colorful = 0;
//      bool color = false;

//      if (Animation.p.Count > 0)
//      {
//        // objects=<int>
//        if (Util.TryParse(Animation.p, "speed", ref speed))
//        {
//          if (speed < 1)
//            speed = 1;
//        }
        

//        if (Util.TryParse(Animation.p, "colorful", ref colorful))
//        {
//          if (colorful!= 1 && colorful!= 0)
//            colorful = 0;
//          if (colorful== 1)
//            color = true;
//          else
//            color = false;
//        }
//        if (Util.TryParse(Animation.p, "Xrotation", ref Xrotation))
//        {
//          if (Xrotation < -1)
//            Xrotation = -1;
//          if (Xrotation > 1)
//            Xrotation = 1;
//        }
//        if (Util.TryParse(Animation.p, "Yrotation", ref Yrotation))
//        {
//          if (Yrotation < -1)
//            Yrotation = -1;
//          if (Yrotation > 1)
//            Yrotation = 1;
//        }
//        if (Util.TryParse(Animation.p, "Zrotation", ref Zrotation))
//        {
//          if (Zrotation < -1)
//            Zrotation = -1;
//          if (Zrotation > 1)
//            Zrotation = 1;
//        }
//      }

//      double t = ((time - start) / (end - start));
      
//      b.Yvalues.Keys.CopyTo(yVals, 0);
//      Array.Sort(yVals);
//      foreach (double yVal in yVals)
//      {
//        Point p = b.Yvalues[yVal];
//        p.RotateVector(AXIS.axisZ, Math.PI * Zrotation * speed * t);
//        p.RotateVector(AXIS.axisY, Math.PI * Yrotation * speed * t);
//        p.RotateVector(AXIS.axisX, Math.PI * Xrotation * speed * t);


//        p.Get2DRepresentation(out double x, out double y, out double radius, c);

//        c.SetColor(b.GetColor(p, color));

//        c.FillDisc((float)x, (float)y, (float)radius);


//        p.RotateVector(AXIS.axisX, -Math.PI * Xrotation * speed * t);
//        p.RotateVector(AXIS.axisY, -Math.PI * Yrotation * speed * t);
//        p.RotateVector(AXIS.axisZ, -Math.PI * Zrotation * speed * t);
//      }
//      //b.MoveShade(speed);

//    }


//    public enum AXIS
//    {
//      axisZ,
//      axisY,
//      axisX
//    }

//    public class ShadePoint : Point
//    {
//      public Point Direction { get; private set; }
//      public double[][] BaseMatrix { get; private set; }

//      public ShadePoint (RandomJames r)
//      {
//        //Direction = new Point();
//        //// sets random direction vector 
//        //Direction.RotateVector(AXIS.axisX, r.RandomDouble(0, Math.PI));
//        //Direction.RotateVector(AXIS.axisY, r.RandomDouble(0, Math.PI));
//        //Direction.Z = 0;
//        //double norm = Math.Sqrt(X*X + Y*Y);
//        //X /= norm;
//        //Y /= norm;

//        //Point vector2 = new Point(X,Y,Z);
//        //vector2.RotateVector(AXIS.axisX, Math.PI / 2);
//        //Point vector3 = new Point(X,Y,Z);
//        //vector3.RotateVector(AXIS.axisY, Math.PI / 2);

//        //BaseMatrix = new double[3][] { new double[3],new double[3], new double[3]};
//        //BaseMatrix[0][0] = X;
//        //BaseMatrix[0][1] = Y;
//        //BaseMatrix[0][2] = Z;
//        //BaseMatrix[1][0] = vector2.X;
//        //BaseMatrix[1][1] = vector2.Y;
//        //BaseMatrix[1][2] = vector2.Z;
//        //BaseMatrix[2][0] = vector3.X;
//        //BaseMatrix[2][1] = vector3.Y;
//        //BaseMatrix[2][2] = vector3.Z;

//      }
//    }

//    [DebuggerDisplay("X:{X}, Y:{Y}, Z:{Z}")]
//    public class Point
//    {
//      public double X { get; set; }
//      public double Y { get; set; }
//      public double Z { get; set; }


//      public Point (double x, double y, double z)
//      {
//        this.X = x;
//        this.Y = y;
//        this.Z = z;
//      }

//      public Point ()
//      {
//        this.X = 1;
//        this.Y = 0;
//        this.Z = 0;
//      }

//      /// <summary>
//      /// Rotates vector that represents point in 3D space
//      /// </summary>
//      /// <param name="axis"> rotation about axisX / axisY / axisZ </param>
//      /// <param name="angle">angle in radians (positive -> anticlockwise / negative -> clockwise </param>
//      public void RotateVector (AXIS axis, double angle)
//      {
//        switch (axis)
//        {
//          case AXIS.axisZ:
//            double tempX = Math.Cos(angle) * this.X - Math.Sin(angle) * this.Y + 0 * this.Z;
//            double tempY = Math.Sin(angle) * this.X + Math.Cos(angle) * this.Y + 0 * this.Z;
//            double tempZ = 0 * this.X + 0 * this.Y + 1 * this.Z;
//            this.X = tempX;
//            this.Y = tempY;
//            this.Z = tempZ;


//            break;
//          case AXIS.axisY:
//            tempX = Math.Cos(angle) * this.X + 0 * this.Y + Math.Sin(angle) * this.Z;
//            tempY = 0 * this.X + 1 * this.Y + 0 * this.Z;
//            tempZ = -Math.Sin(angle) * this.X + 0 * this.Y + Math.Cos(angle) * this.Z;
//            this.X = tempX;
//            this.Y = tempY;
//            this.Z = tempZ;

//            break;
//          case AXIS.axisX:
//            tempX = 1 * this.X + 0 * this.Y + 0 * this.Z;
//            tempY = 0 * this.X + Math.Cos(angle) * this.Y - Math.Sin(angle) * this.Z;
//            tempZ = 0 * this.X + Math.Sin(angle) * this.Y + Math.Cos(angle) * this.Z;
//            this.X = tempX;
//            this.Y = tempY;
//            this.Z = tempZ;

//            break;
//          default:
//            break;
//        }


//      }

//      public void Get2DRepresentation (out double x, out double y, out double radius, Canvas c)
//      {
//        double niceConstantThatRepresentsRatioBetweenHeightAndWidthOfTheImage = ((c.Width/10)*(c.Height/10))/(c.Width+c.Height);
//        double multiplicator = 0.5f;

//        y = (c.Height / 2 + Z * (c.Height / 4));
//        x = (c.Width / 2 + X * (c.Width / 4));
//        radius = (float)(1 + niceConstantThatRepresentsRatioBetweenHeightAndWidthOfTheImage * (Y + 1) * (Y + 1));

//      }
//    }



//    public class Ball
//    {

//      // this is going to be ugly but it will help us later //
//      public Dictionary<double, Point> Yvalues { get; private set; }
//      public double Radius { get; private set; }
//      public List<Point> Points { get; private set; }
//      public ShadePoint shade { get; set; }
//      public float shadowSize;
//      public Ball (int radius, int horizontalLines, int verticalLines, long seed, float shadowSize)
//      {
//        this.shadowSize = shadowSize;
//        Radius = radius;
//        Yvalues = new Dictionary<double, Point>();
//        Points = GeneratePoints(horizontalLines, verticalLines);
        
//        RandomJames r = new RandomJames(seed);
//        shade = new ShadePoint(r);
//        //sets random position on ball //
//        shade.RotateVector(AXIS.axisX, r.RandomDouble(0, Math.PI));
//        shade.RotateVector(AXIS.axisY, r.RandomDouble(0, Math.PI));
//        shade.RotateVector(AXIS.axisZ, r.RandomDouble(0, Math.PI));
//      }

//      List<Point> GeneratePoints (int HL, int VL)
//      {
//        List<Point> points = new List<Point>();
//        for (int i = 1; i < HL; i++)
//        {
//          for (int j = 0; j < VL; j++)
//          {

//            // https://stackoverflow.com/questions/4081898/procedurally-generate-a-sphere-mesh <- credit
//            Point p = new Point(
//              Math.Sin(Math.PI * (double)((double)i/(double)HL)) * Math.Cos(2* Math.PI * (double)((double)j/(double)VL)),
//              Math.Sin(Math.PI * (double)((double)i/(double)HL)) * Math.Sin(2* Math.PI * (double)((double)j/(double)VL)),
//              Math.Cos(Math.PI * (double)((double)i/(double)HL)));

//            points.Add(p);
//            // dont blame me pls
//            while (Yvalues.ContainsKey(p.Y)){
//              p.Y += 0.0000000001;
//            }
//            Yvalues.Add(p.Y, p);
//          }
//        }
//        return points;
//      }
//      //public void MoveShade(double speed)
//      //{
//      //  shade.X = shade.Direction.X*speed + shade.X;
//      //  shade.Y = shade.Direction.Y*speed + shade.Y;
//      //  shade.Z = shade.Direction.Z*speed + shade.Z;
//      //  // norm again
//      //  double norm = Math.Sqrt(shade.X*shade.X + shade.Y*shade.Y + shade.Z*shade.Z);
//      //  shade.X /= norm;
//      //  shade.Y /= norm;
//      //  shade.Z /= norm;
//      //}

//      public Color GetColor (Point p, bool color)
//      {

//        double xDiff = shade.X - p.X;
//        double yDiff = shade.Y - p.Y;
//        double zDiff = shade.Z - p.Z;


//        double distance = Math.Abs(Math.Sqrt(xDiff*xDiff + yDiff*yDiff + zDiff*zDiff));

//        double coeff = 1;
//        int r = 255;
//        int b = 255;
//        int g = 255;
//        if (color)
//        {
//          // could have done something more "nicer", but im too tired to do anything else (3am)
//          r = (int)(Math.Abs((p.X + p.Y + p.Z)/3)*255);
//          g = (int)(Math.Abs((p.X - p.Y + p.Z)/3)*255);
//          b = (int)(Math.Abs((p.X + p.Y - p.Z)/3)*255);
//        }
//        if (distance < shadowSize)
//        {
//          coeff = (double)((distance / shadowSize) * (distance / shadowSize));
//        }

//        return Color.FromArgb(255, (int)(coeff*r), (int)(coeff*g), (int)(coeff*b));
//      }
//    }
//  }
//}
