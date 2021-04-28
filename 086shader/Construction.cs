using System;
using OpenTK;
using Utilities;
using System.Collections.Generic;

namespace Scene3D
{

  public class Curve
  {
    public double xFreq, yFreq, zFreq, thickness, vertexFreq, circEdgeSize;

    public List<Vector3> Vertices { get; private set; }
    public List<List<Vector3>> MeshVertices { get; private set; }
    public int GeneratedFaces { get; private set; }
    public int CurveSegments { get; set; }
    public int MeshCircleCount { get; set; }

    private int CurveCirclesBaseVal;
    public Curve (double xFreq, double yFreq, double zFreq, double thickness, double vertexFreq, double circEdgeSize)
    {
      this.xFreq = xFreq;
      this.yFreq = yFreq;
      this.zFreq = zFreq;
      this.thickness = thickness;
      this.vertexFreq = vertexFreq;
      this.circEdgeSize = circEdgeSize;
      this.CurveCirclesBaseVal = 1000;
      this.CurveSegments = (int)(this.CurveCirclesBaseVal * vertexFreq);
      this.Vertices = new List<Vector3>();
      this.MeshVertices = new List<List<Vector3>>();
      this.MeshCircleCount = (int)(3*(circEdgeSize > 1 ? circEdgeSize : 1));
      this.GeneratedFaces = 1;
    }

    internal void GenerateThingsForScene(Matrix4 m) {
      double fullRotation = 2 * Math.PI;
      double oneSegmentSpacing = fullRotation / CurveSegments;
      double currVal = 0;

      // creates points of the curve itself based on vertexFreq
      for (int i = 0; i < CurveSegments; i++) {
        var cVector = new Vector3(
          (float)Math.Sin(currVal * xFreq),
          (float)Math.Sin(currVal * yFreq),
          (float)Math.Sin(currVal * zFreq));
        Vertices.Add(cVector);
        currVal += oneSegmentSpacing;
      }

      CreateVerticesAroundPoints();
      return;
    }

    private void CreateVerticesAroundPoints() {
      double oneSegmentSpacing = (2 * Math.PI) / MeshCircleCount;
      for (int i = 0; i < Vertices.Count;i++) {
        List<Vector3> cV = new List<Vector3>();
        double currVal = 0;
        for (int j = 0; j < MeshCircleCount; j++) {
          cV.Add(GetRotatedVectorThatRepresentsPartOfCircle(
            Vertices[(i == 0) ? Vertices.Count-1 : i-1],
            Vertices[i],
            Vertices[(i == Vertices.Count-1) ? 0 : i+1],
            currVal));
          currVal += oneSegmentSpacing;
        }
        MeshVertices.Add(cV);
      }
    }

    internal Vector3 GetRotatedVectorThatRepresentsPartOfCircle(Vector3 prev, Vector3 curr, Vector3 next, double currVal) {
      //finding line (represented as vector) such that both prev and next vectors point at it
      //Vector3 line = new Vector3(next.X - prev.X, next.Y - prev.Y, next.Z - prev.Z);
      /*// i have two(prev and next (actually three ^)) vectors in same plane, now i have to
      // gramschmidt (map) current to that plane
      Vector3 u1 = prev;
      float fractionu1next = (u1.X * next.X + u1.Y * next.Y + u1.Z * next.Z) /
                        (u1.X * u1.X + u1.Y * u1.Y + u1.Z * u1.Z);
      Vector3 u2 = next - (u1.scale(fractionu1next,fractionu1next,fractionu1next));
      Vector3 baseVector1 = u1.Normalized();
      Vector3 baseVector2 = u2.Normalized();
      // now we have two base vectors that are perpendicular, we can now derive and project
      // curr vector to the plane created by these two vectors
      float fractionu1curr = (u1.X * curr.X + u1.Y * curr.Y + u1.Z * curr.Z) /
                            (u1.X * u1.X + u1.Y * u1.Y + u1.Z * u1.Z);
      Vector3 u1curr = u1.scale(fractionu1curr, fractionu1curr, fractionu1curr);
      float fractionu2curr = (u2.X * curr.X + u2.Y * curr.Y + u2.Z * curr.Z) /
                            (u2.X * u2.X + u2.Y * u2.Y + u2.Z * u2.Z);
      Vector3 u2curr = u2.scale(fractionu2curr, fractionu2curr, fractionu2curr);
      Vector3 projectedCurr = u1curr + u2curr;

      // after this i realised that i only need angle between prev and next, that gives me the tilt
      // to get the angle i just need vector between prev and next, using pythagorean gives me the angle 
      */
      //float sizeSquared = line.X * line.X + line.Y * line.Y + line.Z * line.Z;
      //double xAngle = Math.Asin(curr.X);
      //double yAngle = Math.Asin(curr.Z);
      //double zAngle = Math.Asin(curr.Y);

      //// i tried...
      //return res.rotate(xAngle,'y').rotate(yAngle,'z').rotate(zAngle,'x');
      Matrix3 baseMatrix;
      if (curr.Length != 0) {
        Vector3 v2 = curr.rotate(Math.PI / 2, 'x');
        Vector3 v3 = curr.rotate(Math.PI / 2, 'y');
        baseMatrix = new Matrix3(v2.Normalized(), curr.Normalized(), v3.Normalized());
      }
      else {
        baseMatrix = new Matrix3(new Vector3(1,0,0), new Vector3(0,1,0), new Vector3(0,0,1));
      }
      // curr, v2 and v3 make orthonormal base, x/z are v2 and v3, y(height) is original curr
      Vector3 v = new Vector3((float)(thickness*Math.Sin(currVal)), 0,(float)(thickness*Math.Cos(currVal) ));

      Vector3 result = new Vector3(
        curr.X + (baseMatrix.M11 * v.X + baseMatrix.M12 * v.Y + baseMatrix.M13 * v.Z),
        curr.Y + (baseMatrix.M21 * v.X + baseMatrix.M22 * v.Y + baseMatrix.M23 * v.Z),
        curr.Z + (baseMatrix.M31 * v.X + baseMatrix.M32 * v.Y + baseMatrix.M33 * v.Z));

      //willHelpUsLater.Normalize();
      return result;
    }
    internal void RenderSceneJustCurve(SceneBrep scene) {
      for (int i = 0; i < Vertices.Count; i++) {
        scene.AddVertex(Vertices[i]);
        if (i > 0) scene.AddLine(i - 1, i);
      }
      scene.AddLine(0, Vertices.Count - 1);
    }

    internal void RenderSceneSimpleTube(SceneBrep scene) {
      int total = -1;
      for(int i = 0; i < MeshVertices.Count; i++) {
        for (int j = 0; j < MeshVertices[i].Count; j++) {  
          scene.AddVertex(MeshVertices[i][j]);
          total++;
        }
        //AddToScene(scene, MeshVertices[i], total);
        if (i > 0) {
          ConnectTwoCircles(scene, MeshVertices, total, i);
        }
      }
      ConnectTwoCircles(scene, MeshVertices, total, 0);
    }

    void ConnectTwoCircles(SceneBrep scene, List<List<Vector3>> MeshV, int total, int index) {
      List<int> firstCirc = new List<int>();
      List<int> secCirc = new List<int>();
      int count = MeshV[index].Count;
      if (index != 0) {
        for (int i = total; i > total - count; i--) {
          firstCirc.Add(i);
          secCirc.Add(i - count);
        }
        //for (int i = 0; i < firstCirc.Count; i++) {
        //  scene.AddLine(firstCirc[i], secCirc[i]);
        //  scene.AddLine(firstCirc[i], secCirc[(i + 1 == secCirc .Count ? 0 : i+1)]);
        //}
      }
      else {
        for (int i = total; i > total - count; i--) {
          firstCirc.Add(i);
          secCirc.Add(total-i+count);
        }
        //for (int i = 0; i < firstCirc.Count; i++) {
        //  scene.AddLine(firstCirc[i], secCirc[i]);
        //  scene.AddLine(firstCirc[i], secCirc[(i + 1 == secCirc.Count ? 0 : i+1)]);
        //}
      }

      for (int i = 0; i < firstCirc.Count; i++) {
        scene.AddTriangle(firstCirc[i], secCirc[i], secCirc[(i + 1 == secCirc.Count ? 0 : i + 1)]);
        this.GeneratedFaces++;
        scene.AddTriangle(
          firstCirc[(i + 1 == firstCirc.Count ? 0 : i + 1)],
          firstCirc[i],
          secCirc[(i + 1 == secCirc.Count ? 0 : i + 1)]);
        this.GeneratedFaces++;
      }
    }

    void AddToScene(SceneBrep scene, List<Vector3> v, int total) {
      for(int i = total; i > total-v.Count+1; i--) {
        scene.AddLine(i, i - 1);
      }
    }
  }

  public class Construction
  {
    #region Form initialization

    /// <summary>
    /// Optional data initialization.
    /// </summary>
    /// 
    public static double xFreq, yFreq, zFreq;
    public static double thickness, vertexFreq, circEdgeSize;
    public static void InitParams (out string name, out string param, out string tooltip)
    {
      name    = "Lukáš Salak";
      param   = "xFreq = 2, yFreq = 5, zFreq = 3, thickness = 0.08, vertexFreq = 0.5, circEdgeSize = 4";
      tooltip = "";
        }

    #endregion

    #region Instance data

    // Instance data.
    

    #endregion

    #region Construction

    public Construction ()
    {
      // Any one-time initialization code goes here..
    }

    private void parseParams (string param)
    {

      Dictionary<string, string> p = Util.ParseKeyValueList( param );
      if (p.Count > 0) {
        Util.TryParse(p, "xFreq", ref xFreq);
        Util.TryParse(p, "yFreq", ref yFreq);
        Util.TryParse(p, "zFreq", ref zFreq);
        Util.TryParse(p, "thickness", ref thickness);
        Util.TryParse(p, "vertexFreq", ref vertexFreq);
        Util.TryParse(p, "circEdgeSize", ref circEdgeSize);
      }
    }
    #endregion

    #region Mesh construction

    /// <summary>
    /// Construct a new Brep solid (preferebaly closed = regular one).
    /// </summary>
    /// <param name="scene">B-rep scene to be modified</param>
    /// <param name="m">Transform matrix (object-space to world-space)</param>
    /// <param name="param">Shape parameters if needed</param>
    /// <returns>Number of generated faces (0 in case of failure)</returns>
    public int AddMesh (SceneBrep scene, Matrix4 m, string param)
    {
      parseParams(param);
      scene.Reserve(10000000);
      Curve c = new Curve(xFreq,yFreq, zFreq,thickness,vertexFreq,circEdgeSize);
      c.GenerateThingsForScene(m);
      c.RenderSceneSimpleTube(scene);
     // c.RenderSceneJustCurve(scene);

      return c.GeneratedFaces;
    }


    

    #endregion
  }
  public static class Transforms {
    public static Vector3[] scale(Vector3[] to_scale, float x, float y, float z) {
      Vector3[] scaled = new Vector3[to_scale.Length];
      MyMatrix4 mat = MyMatrix4.getScaleMatrix(x, y, z);
      for (int i = 0; i < to_scale.Length; i++)
        scaled[i] = mat * to_scale[i];

      return scaled;
    }

    public static Vector3 scale(this Vector3 to_scale, float x, float y, float z) {
      MyMatrix4 mat = MyMatrix4.getScaleMatrix(x, y, z);
      return mat * to_scale;
    }

    public static Vector3[] translate(Vector3[] to_translate, float x, float y, float z) {
      Vector3[] translated = new Vector3[to_translate.Length];
      MyMatrix4 mat = MyMatrix4.getTranslationMatrix(x, y, z);
      for (int i = 0; i < to_translate.Length; i++)
        translated[i] = mat * to_translate[i];

      return translated;
    }

    public static Vector3 translate(Vector3 to_translate, float x, float y, float z) {
      MyMatrix4 mat = MyMatrix4.getTranslationMatrix(x, y, z);
      return mat * to_translate;
    }

    public static Vector3[] rotate(Vector3[] to_rotate, double angle, char axis) {
      Vector3[] rotated = new Vector3[to_rotate.Length];
      MyMatrix4 mat = MyMatrix4.getSingleRotationMatrix(angle, axis);
      for (int i = 0; i < to_rotate.Length; i++)
        rotated[i] = mat * to_rotate[i];

      return rotated;
    }

    public static Vector3 rotate(this Vector3 to_rotate, double angle, char axis) {
      MyMatrix4 mat = MyMatrix4.getSingleRotationMatrix(angle, axis);
      return mat * to_rotate;
    }

    public static Vector3[] flip(Vector3[] to_flip, char axis) {
      Vector3[] flipped = new Vector3[to_flip.Length];
      for (int i = 0; i < to_flip.Length; i++) {
        flipped[i] = to_flip[i];
        switch (axis) {
          case 'x':
            flipped[i].X = -flipped[i].X;
            break;
          case 'y':
            flipped[i].Y = -flipped[i].Y;
            break;
          case 'z':
            flipped[i].Z = -flipped[i].Z;
            break;
        }
      }

      return flipped;
    }

    public static Vector3[] reverse(Vector3[] to_reverse) {
      Vector3[] reversed = new Vector3[to_reverse.Length];
      for (int i = 0; i < to_reverse.Length; i++)
        reversed[i] = to_reverse[to_reverse.Length - 1 - i];

      return reversed;
    }
  }

  public class MyMatrix4 {
    public float[,] Lattice { get; set; }

    public MyMatrix4() {
      Lattice = new float[4, 4];
    }

    public static MyMatrix4 getTranslationMatrix(float x, float y, float z) {
      MyMatrix4 m = getOneMatrix();
      m.Lattice[0, 3] = x;
      m.Lattice[1, 3] = y;
      m.Lattice[2, 3] = z;

      return m;
    }

    public static MyMatrix4 getScaleMatrix(float x, float y, float z) {
      MyMatrix4 m = getOneMatrix();
      m.Lattice[0, 0] = x;
      m.Lattice[1, 1] = y;
      m.Lattice[2, 2] = z;

      return m;
    }

    public static MyMatrix4 getSingleRotationMatrix(double angle, char axis) {
      MyMatrix4 m = getOneMatrix();
      switch (axis) {
        case 'x':
          m.Lattice[1, 1] = (float)Math.Cos(angle);
          m.Lattice[1, 2] = (float)-Math.Sin(angle);
          m.Lattice[2, 1] = (float)Math.Sin(angle);
          m.Lattice[2, 2] = (float)Math.Cos(angle);
          break;
        case 'y':
          m.Lattice[0, 0] = (float)Math.Cos(angle);
          m.Lattice[0, 2] = (float)Math.Sin(angle);
          m.Lattice[2, 0] = (float)-Math.Sin(angle);
          m.Lattice[2, 2] = (float)Math.Cos(angle);
          break;
        case 'z':
          m.Lattice[0, 0] = (float)Math.Cos(angle);
          m.Lattice[0, 1] = (float)-Math.Sin(angle);
          m.Lattice[1, 0] = (float)Math.Sin(angle);
          m.Lattice[1, 1] = (float)Math.Cos(angle);
          break;
      }

      return m;
    }

    public static MyMatrix4 getOneMatrix() {
      MyMatrix4 m = new MyMatrix4();
      for (int i = 0; i < 4; i++)
        m.Lattice[i, i] = 1;

      return m;
    }

    public static Vector3 operator *(MyMatrix4 m, Vector3 v) {
      float a = v.X * m.Lattice[0, 0] + v.Y * m.Lattice[0, 1] + v.Z * m.Lattice[0, 2] + m.Lattice[0, 3];
      float b = v.X * m.Lattice[1, 0] + v.Y * m.Lattice[1, 1] + v.Z * m.Lattice[1, 2] + m.Lattice[1, 3];
      float c = v.X * m.Lattice[2, 0] + v.Y * m.Lattice[2, 1] + v.Z * m.Lattice[2, 2] + m.Lattice[2, 3];
      float d = v.X * m.Lattice[3, 0] + v.Y * m.Lattice[3, 1] + v.Z * m.Lattice[3, 2] + m.Lattice[3, 3];

      return new Vector3(a, b, c);
    }
  }
}

// using System;
// using OpenTK;
// using Utilities;
// using System.Collections.Generic;

// namespace Scene3D
// {
  // public class Construction
  // {
    // #region Form initialization

    // /// <summary>
    // /// Optional data initialization.
    // /// </summary>
    // public static void InitParams (out string name, out string param, out string tooltip)
    // {
      // name    = "Matúš Goliaš";
      // param   = "level=5, size=0.2, cdetail=16, radius=0.75, bdetail=9";
      // tooltip = "level=recursion-level, size=relative pipe cross section\ncdetail=cross section divisions, radius=relative knee-radius\nbdetail=knee radius";
    // }

    // #endregion

    // #region Instance data

    // // Instance data.
    // private Curve curve;
    // private float radius;
    // private int cdetail;
    // private int bdetail;
    // private int level;
    // private float size;
    // private Vector3[] last_circle;
    // private int[] connect_indices;

    // #endregion

    // #region Construction

    // public Construction ()
    // {
      // // Any one-time initialization code goes here..
    // }

    // private void parseParams (string param)
    // {
      // level = 2;
      // radius = 0.25f;
      // size = 0.25f;
      // cdetail = 5;
      // bdetail = 5;
      // Dictionary<string, string> p = Util.ParseKeyValueList( param );
      // if (p.Count > 0)
      // {
        // if (Util.TryParse(p, "level", ref level) && level < 1)
        // {
          // level = 1;
        // }
        // if (Util.TryParse(p, "size", ref size) && (size < 0 || size > 0.49f))
        // {
          // size = 0.25f;
        // }
        // if (Util.TryParse(p, "radius", ref radius))
        // {
          // if (radius < 0)
            // radius = 0.25f;
          // else
          // {
            // if (radius * 2 > size)
            // {
              // radius = 2 * size - 0.05f;
              // if (radius < 0)
                // radius = 0.05f;
            // }
          // }
        // }
        // if (Util.TryParse(p, "cdetail", ref cdetail) && cdetail < 3)
        // {
          // cdetail = 3;
        // }
        // if (Util.TryParse(p, "bdetail", ref bdetail) && bdetail < 3)
        // {
          // bdetail = 3;
        // }
      // }
    // }

    // #endregion

    // #region Mesh construction

    // /// <summary>
    // /// Construct a new Brep solid (preferebaly closed = regular one).
    // /// </summary>
    // /// <param name="scene">B-rep scene to be modified</param>
    // /// <param name="m">Transform matrix (object-space to world-space)</param>
    // /// <param name="param">Shape parameters if needed</param>
    // /// <returns>Number of generated faces (0 in case of failure)</returns>
    // public int AddMesh (SceneBrep scene, Matrix4 m, string param)
    // {
      // parseParams(param);
      // scene.Reserve(10000000);

      // curve = Curve.createCurve(level);

      // firstCircle(curve.current_vertices[0], curve.current_vertices[1], scene, m);
      // int count = 0;

      // for (int i = 1; i < curve.current_vertices.Length - 1; i++)
      // {
        // if (isBent(curve.current_vertices[i - 1], curve.current_vertices[i], curve.current_vertices[i + 1]))
          // count += createBentPipe(scene, m, curve.current_vertices[i - 1], curve.current_vertices[i], curve.current_vertices[i + 1]);
        // else
        // {
          // Vector3 dif = curve.current_vertices[i] - curve.current_vertices[i - 1];
          // for (int j = 0; j < cdetail; j++)
            // last_circle[j] = last_circle[j] + dif;
        // }
      // }

      // int l = curve.current_vertices.Length;
      // count += lastCircle(curve.current_vertices[l - 2], curve.current_vertices[l - 1], scene, m);

      // return count;
    // }

    // private bool isBent (Vector3 v1, Vector3 v2, Vector3 v3)
    // {
      // Vector3 a = v2 - v1;
      // Vector3 b = v3 - v2;
      // a.X = Math.Abs(a.X);
      // a.Y = Math.Abs(a.Y);
      // a.Z = Math.Abs(a.Z);
      // b.X = Math.Abs(b.X);
      // b.Y = Math.Abs(b.Y);
      // b.Z = Math.Abs(b.Z);
      // int al, bl;

      // if (a.Y > a.X)
        // al = (a.Z > a.Y) ? 3 : 2;
      // else
        // al = (a.Z > a.X) ? 3 : 1;

      // if (b.Y > b.X)
        // bl = (b.Z > b.Y) ? 3 : 2;
      // else
        // bl = (b.Z > b.X) ? 3 : 1;

      // return al != bl;
    // }

    // private void firstCircle (Vector3 v1, Vector3 v2, SceneBrep scene, Matrix4 m)
    // {
      // Vector3 dif = (v2 - v1) * size;
      // Vector3 rdif = (v2 - v1) * (1 - size);
      // last_circle = createCircleAround(v1, v2);
      // float s = radius / curve.current_count;
      // last_circle = Transforms.scale(last_circle, s, s, s);
      // connect_indices = new int[cdetail];
      // Vector3 center = getCenter(last_circle);

      // for (int i = 0; i < cdetail; i++)
      // {
        // connect_indices[i] = scene.AddVertex(Vector3.TransformPosition(last_circle[i] + v1, m));
        // scene.SetNormal(connect_indices[i], Vector3.TransformVector(last_circle[i] - center, m).Normalized());
        // float r = (last_circle[ i ].X + v1.X + 1.1f) / 2.5f;
        // float g = (last_circle[ i ].Y + v1.Y + 1.1f) / 2.5f;
        // float b = (last_circle[ i ].Z + v1.Z + 1.1f) / 2.5f;
        // scene.SetColor(connect_indices[i], new Vector3(r, g, b));
      // }

      // for (int i = 0; i < cdetail; i++)
        // last_circle[i] += v1 + rdif;
    // }

    // private int lastCircle (Vector3 v1, Vector3 v2, SceneBrep scene, Matrix4 m)
    // {
      // _ = createCircleAround(v1, v2);
      // Vector3 center = getCenter(last_circle);
      // Vector3 dif = (v2 - v1) * size;
      // _ = (v2 - v1) * (1 - size);
      // int[] v = new int[cdetail];

      // for (int i = 0; i < cdetail; i++)
      // {
        // v[i] = scene.AddVertex(Vector3.TransformPosition(last_circle[i] + dif, m));
        // scene.SetNormal(v[i], Vector3.TransformVector(last_circle[i] - center, m).Normalized());
        // float r = (last_circle[ i ].X + 1.1f) / 2.5f;
        // float g = (last_circle[ i ].Y + 1.1f) / 2.5f;
        // float b = (last_circle[ i ].Z + 1.1f) / 2.5f;
        // scene.SetColor(v[i], new Vector3(r, g, b));
      // }

      // for (int i = 0; i < cdetail - 1; i++)
      // {
        // scene.AddTriangle(connect_indices[i], v[i], connect_indices[i + 1]);
        // scene.AddTriangle(connect_indices[i + 1], v[i], v[i + 1]);
      // }

      // scene.AddTriangle(connect_indices[cdetail - 1], v[cdetail - 1], connect_indices[0]);
      // scene.AddTriangle(connect_indices[0], v[cdetail - 1], v[0]);

      // return cdetail * 2;
    // }

    // private Vector3[] createCircleAround (Vector3 v1, Vector3 v2)
    // {
      // Vector3 rot_point = v2 - v1;
      // rot_point.X = Math.Abs(rot_point.X);
      // rot_point.Y = Math.Abs(rot_point.Y);
      // rot_point.Z = Math.Abs(rot_point.Z);

      // if (rot_point.Y > rot_point.X)
      // {
        // if (rot_point.Z > rot_point.Y)
          // return createCircleXY(cdetail);
        // return createCircleXZ(cdetail);
      // }

      // if (rot_point.Z > rot_point.X)
        // return createCircleXY(cdetail);
      // return createCircleYZ(cdetail);
    // }

    // private Vector3[] createCircleXY (int detail)
    // {
      // Vector3[] field = new Vector3[detail];
      // double delta_angle = Math.PI * 2.0 / detail;
      // field[0] = new Vector3(0, 1, 0);
      // for (int i = 1; i < detail; i++)
        // field[i] = Transforms.rotate(field[0], delta_angle * i, 'z');

      // return field;
    // }

    // private Vector3[] createCircleXZ (int detail)
    // {
      // Vector3[] field = new Vector3[detail];
      // double delta_angle = Math.PI * 2.0 / detail;
      // field[0] = new Vector3(0, 0, 1);
      // for (int i = 1; i < detail; i++)
        // field[i] = Transforms.rotate(field[0], delta_angle * i, 'y');

      // return field;
    // }

    // private Vector3[] createCircleYZ (int detail)
    // {
      // Vector3[] field = new Vector3[detail];
      // double delta_angle = Math.PI * 2.0 / detail;
      // field[0] = new Vector3(0, 0, 1);
      // for (int i = 1; i < detail; i++)
        // field[i] = Transforms.rotate(field[0], delta_angle * i, 'x');

      // return field;
    // }

    // private float longestSide (Vector3 v)
    // {
      // float x = Math.Abs(v.X);
      // float y = Math.Abs(v.Y);
      // float z = Math.Abs(v.Z);
      // return Math.Max(x, Math.Max(y, z));
    // }

    // private int createBentPipe (SceneBrep scene, Matrix4 m, Vector3 v1, Vector3 v2, Vector3 v3)
    // {
      // Vector3[] centers = new Vector3[bdetail];
      // Vector3 oned = (v1 - v2) * size;
      // Vector3 twod = (v3 - v2) * size;
      // Vector3 rtwod = (v3 - v2) * (1 - size);
      // Vector3 rot_point = v2 + oned + twod;
      // Vector3 test_point = oned + twod;
      // test_point.X = Math.Abs(test_point.X);
      // test_point.Y = Math.Abs(test_point.Y);
      // test_point.Z = Math.Abs(test_point.Z);
      // double angle = -(Math.PI / 2) / (bdetail - 1);
      // char axis;

      // if (test_point.Y < test_point.X)
        // axis = (test_point.Z < test_point.Y) ? 'z' : 'y';
      // else
        // axis = (test_point.Z < test_point.X) ? 'z' : 'x';

      // MyMatrix4 mat = MyMatrix4.getSingleRotationMatrix(angle, axis);
      // float l = longestSide(last_circle[0] - v1);
      // Vector3 sp = rotateAroundPoint(last_circle[0], rot_point, mat);
      // float l2 = longestSide(sp - v1);
      // if (l2 < l)
        // mat = MyMatrix4.getSingleRotationMatrix(-angle, axis);
      // Vector3[] field = new Vector3[bdetail * cdetail];
      // Array.Copy(last_circle, field, cdetail);
      // centers[0] = getCenter(last_circle);
      // Vector3[] new_circle = new Vector3[cdetail];

      // for (int i = 1; i < bdetail; i++)
      // {
        // Array.ConstrainedCopy(field, (i - 1) * cdetail, new_circle, 0, cdetail);
        // new_circle = rotateAroundPoint(new_circle, rot_point, mat);
        // Array.ConstrainedCopy(new_circle, 0, field, i * cdetail, cdetail);
        // centers[i] = getCenter(new_circle);
      // }

      // Array.ConstrainedCopy(field, (bdetail - 1) * cdetail, last_circle, 0, cdetail);
      // for (int i = 0; i < cdetail; i++)
        // last_circle[i] += rtwod - twod;
      // int[] v = new int[field.Length];

      // for (int i = 0; i < field.Length; i++)
      // {
        // v[i] = scene.AddVertex(Vector3.TransformPosition(field[i], m));
        // scene.SetNormal(v[i], Vector3.TransformVector(field[i] - centers[i / cdetail], m).Normalized());
        // float r = (field[i].X + 1.1f) / 2.5f;
        // float g = (field[i].Y + 1.1f) / 2.5f;
        // float b = (field[i].Z + 1.1f) / 2.5f;
        // scene.SetColor(v[i], new Vector3(r, g, b));
      // }

      // for (int i = 0; i < cdetail - 1; i++)
      // {
        // scene.AddTriangle(connect_indices[i], v[i], connect_indices[i + 1]);
        // scene.AddTriangle(connect_indices[i + 1], v[i], v[i + 1]);
      // }

      // scene.AddTriangle(connect_indices[cdetail - 1], v[cdetail - 1], connect_indices[0]);
      // scene.AddTriangle(connect_indices[0], v[cdetail - 1], v[0]);

      // for (int i = 0; i < (bdetail - 1); i++)
      // {
        // for (int j = 0; j < cdetail - 1; j++)
        // {
          // int a = i * cdetail + j;
          // scene.AddTriangle(v[a], v[a + cdetail], v[a + 1]);
          // scene.AddTriangle(v[a + 1], v[a + cdetail], v[a + cdetail + 1]);
        // }
        // int b = i * cdetail + cdetail - 1;
        // scene.AddTriangle(v[b], v[b + cdetail], v[b - cdetail + 1]);
        // scene.AddTriangle(v[b - cdetail + 1], v[b + cdetail], v[b + 1]);
      // }
      // Array.ConstrainedCopy(v, (bdetail - 1) * cdetail, connect_indices, 0, cdetail);

      // return bdetail * cdetail * 2 + cdetail * 2;
    // }

    // private Vector3 getCenter (Vector3[] circle)
    // {
      // Vector3 v = new Vector3();
      // foreach (Vector3 c in circle)
        // v += c;

      // v = v / circle.Length;

      // return v;
    // }

    // private Vector3[] rotateAroundPoint (Vector3[] to_rotate, Vector3 point, MyMatrix4 mat)
    // {
      // Vector3[] field = new Vector3[to_rotate.Length];
      // for (int i = 0; i < to_rotate.Length; i++)
        // field[i] = mat * (to_rotate[i] - point) + point;

      // return field;
    // }

    // private Vector3 rotateAroundPoint (Vector3 to_rotate, Vector3 point, MyMatrix4 mat)
    // {
      // return mat * (to_rotate - point) + point;
    // }

    // static class Transforms
    // {
      // public static Vector3[] scale (Vector3[] to_scale, float x, float y, float z)
      // {
        // Vector3[] scaled = new Vector3[to_scale.Length];
        // MyMatrix4 mat = MyMatrix4.getScaleMatrix(x, y, z);
        // for (int i = 0; i < to_scale.Length; i++)
          // scaled[i] = mat * to_scale[i];

        // return scaled;
      // }

      // public static Vector3 scale (Vector3 to_scale, float x, float y, float z)
      // {
        // MyMatrix4 mat = MyMatrix4.getScaleMatrix(x, y, z);
        // return mat * to_scale;
      // }

      // public static Vector3[] translate (Vector3[] to_translate, float x, float y, float z)
      // {
        // Vector3[] translated = new Vector3[to_translate.Length];
        // MyMatrix4 mat = MyMatrix4.getTranslationMatrix(x, y, z);
        // for (int i = 0; i < to_translate.Length; i++)
          // translated[i] = mat * to_translate[i];

        // return translated;
      // }

      // public static Vector3 translate (Vector3 to_translate, float x, float y, float z)
      // {
        // MyMatrix4 mat = MyMatrix4.getTranslationMatrix(x, y, z);
        // return mat * to_translate;
      // }

      // public static Vector3[] rotate (Vector3[] to_rotate, double angle, char axis)
      // {
        // Vector3[] rotated = new Vector3[to_rotate.Length];
        // MyMatrix4 mat = MyMatrix4.getSingleRotationMatrix(angle, axis);
        // for (int i = 0; i < to_rotate.Length; i++)
          // rotated[i] = mat * to_rotate[i];

        // return rotated;
      // }

      // public static Vector3 rotate (Vector3 to_rotate, double angle, char axis)
      // {
        // MyMatrix4 mat = MyMatrix4.getSingleRotationMatrix(angle, axis);
        // return mat * to_rotate;
      // }

      // public static Vector3[] flip (Vector3[] to_flip, char axis)
      // {
        // Vector3[] flipped = new Vector3[to_flip.Length];
        // for (int i = 0; i < to_flip.Length; i++)
        // {
          // flipped[i] = to_flip[i];
          // switch (axis)
          // {
            // case 'x':
              // flipped[i].X = -flipped[i].X;
              // break;
            // case 'y':
              // flipped[i].Y = -flipped[i].Y;
              // break;
            // case 'z':
              // flipped[i].Z = -flipped[i].Z;
              // break;
          // }
        // }

        // return flipped;
      // }

      // public static Vector3[] reverse (Vector3[] to_reverse)
      // {
        // Vector3[] reversed = new Vector3[ to_reverse.Length ];
        // for (int i = 0; i < to_reverse.Length; i++)
          // reversed[i] = to_reverse[to_reverse.Length - 1 - i];

        // return reversed;
      // }
    // }

    // class Curve
    // {
      // public int current_level = 1;
      // public int current_count = 1;
      // public Vector3[] current_vertices;

      // public static Curve createCurve (int level)
      // {
        // Curve curve = new Curve();
        // curve.firstLevel();
        // for (int i = 1; i < level; i++)
        // {
          // curve.current_level++;
          // curve.nextLevel();
        // }

        // return curve;
      // }

      // private void firstLevel ()
      // {
        // Vector3[] field = new Vector3[8];
        // field[0] = new Vector3(-1, -1, -1);
        // field[1] = new Vector3(-1, 1, -1);
        // field[2] = new Vector3(-1, 1, 1);
        // field[3] = new Vector3(-1, -1, 1);
        // field[4] = new Vector3(1, -1, 1);
        // field[5] = new Vector3(1, 1, 1);
        // field[6] = new Vector3(1, 1, -1);
        // field[7] = new Vector3(1, -1, -1);
        // field = Transforms.rotate(field, Math.PI / 2, 'y');
        // field = Transforms.rotate(field, -Math.PI / 2, 'x');
        // field = Transforms.reverse(field);

        // current_vertices = field;
      // }

      // private void nextLevel ()
      // {
        // Vector3[] field = new Vector3[current_vertices.Length * 8];

        // float sz = (float)(current_count) / (current_count * 2 + 1);
        // float dist = (float)(current_count + 1) / (current_count * 2 + 1);
        // current_count = current_count * 2 + 1;

        // // Segment 1 of 8.
        // Vector3[] v1;
        // v1 = Transforms.scale(current_vertices, sz, sz, sz);
        // v1 = Transforms.rotate(v1, Math.PI / 2, 'x');
        // v1 = Transforms.reverse(v1);
        // v1 = Transforms.translate(v1, -dist, -dist, dist);

        // // Segment 2 of 8.
        // Vector3[] v2;
        // v2 = Transforms.scale(current_vertices, sz, sz, sz);
        // v2 = Transforms.rotate(v2, Math.PI / 2, 'z');
        // v2 = Transforms.reverse(v2);
        // v2 = Transforms.translate(v2, -dist, -dist, -dist);

        // // Segment 3 of 8.
        // Vector3[] v3;
        // v3 = Transforms.scale(current_vertices, sz, sz, sz);
        // v3 = Transforms.translate(v3, dist, -dist, -dist);

        // // Segment 4 of 8.
        // Vector3[] v4;
        // v4 = Transforms.scale(current_vertices, sz, sz, sz);
        // v4 = Transforms.rotate(v4, -Math.PI / 2, 'z');
        // v4 = Transforms.rotate(v4, -Math.PI / 2, 'x');
        // v4 = Transforms.translate(v4, dist, -dist, dist);

        // // Segment 5 of 8.
        // Vector3[] v5;
        // v5 = Transforms.flip(v4, 'y');

        // // Segment 6 of 8.
        // Vector3[] v6;
        // v6 = Transforms.flip(v3, 'y');

        // // Segment 7 of 8.
        // Vector3[] v7;
        // v7 = Transforms.flip(v2, 'y');

        // // Segment 8 of 8.
        // Vector3[] v8;
        // v8 = Transforms.flip(v1, 'y');
        // for (int i = 0; i < current_vertices.Length; i++)
        // {
          // int pl = current_vertices.Length;
          // field[i] = v1[i];
          // field[pl + i] = v2[i];
          // field[pl * 2 + i] = v3[i];
          // field[pl * 3 + i] = v4[i];
          // field[pl * 4 + i] = v5[pl - 1 - i];
          // field[pl * 5 + i] = v6[pl - 1 - i];
          // field[pl * 6 + i] = v7[pl - 1 - i];
          // field[pl * 7 + i] = v8[pl - 1 - i];
        // }

        // current_vertices = field;
      // }
    // }

    // class MyMatrix4
    // {
      // public float[,] Lattice { get; set; }

      // public MyMatrix4 ()
      // {
        // Lattice = new float[4, 4];
      // }

      // public static MyMatrix4 getTranslationMatrix (float x, float y, float z)
      // {
        // MyMatrix4 m = getOneMatrix();
        // m.Lattice[0, 3] = x;
        // m.Lattice[1, 3] = y;
        // m.Lattice[2, 3] = z;

        // return m;
      // }

      // public static MyMatrix4 getScaleMatrix (float x, float y, float z)
      // {
        // MyMatrix4 m = getOneMatrix();
        // m.Lattice[0, 0] = x;
        // m.Lattice[1, 1] = y;
        // m.Lattice[2, 2] = z;

        // return m;
      // }

      // public static MyMatrix4 getSingleRotationMatrix (double angle, char axis)
      // {
        // MyMatrix4 m = getOneMatrix();
        // switch (axis)
        // {
          // case 'x':
            // m.Lattice[1, 1] = (float)Math.Cos(angle);
            // m.Lattice[1, 2] = (float)-Math.Sin(angle);
            // m.Lattice[2, 1] = (float)Math.Sin(angle);
            // m.Lattice[2, 2] = (float)Math.Cos(angle);
            // break;
          // case 'y':
            // m.Lattice[0, 0] = (float)Math.Cos(angle);
            // m.Lattice[0, 2] = (float)Math.Sin(angle);
            // m.Lattice[2, 0] = (float)-Math.Sin(angle);
            // m.Lattice[2, 2] = (float)Math.Cos(angle);
            // break;
          // case 'z':
            // m.Lattice[0, 0] = (float)Math.Cos(angle);
            // m.Lattice[0, 1] = (float)-Math.Sin(angle);
            // m.Lattice[1, 0] = (float)Math.Sin(angle);
            // m.Lattice[1, 1] = (float)Math.Cos(angle);
            // break;
        // }

        // return m;
      // }

      // public static MyMatrix4 getOneMatrix ()
      // {
        // MyMatrix4 m = new MyMatrix4();
        // for (int i = 0; i < 4; i++)
          // m.Lattice[i, i] = 1;

        // return m;
      // }

      // public static Vector3 operator * (MyMatrix4 m, Vector3 v)
      // {
        // float a = v.X * m.Lattice[0, 0] + v.Y * m.Lattice[0, 1] + v.Z * m.Lattice[0, 2] + m.Lattice[0, 3];
        // float b = v.X * m.Lattice[1, 0] + v.Y * m.Lattice[1, 1] + v.Z * m.Lattice[1, 2] + m.Lattice[1, 3];
        // float c = v.X * m.Lattice[2, 0] + v.Y * m.Lattice[2, 1] + v.Z * m.Lattice[2, 2] + m.Lattice[2, 3];
        // float d = v.X * m.Lattice[3, 0] + v.Y * m.Lattice[3, 1] + v.Z * m.Lattice[3, 2] + m.Lattice[3, 3];

        // return new Vector3(a, b, c);
      // }
    // }

    // #endregion
  // }
// }
