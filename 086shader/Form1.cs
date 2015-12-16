using System;
using System.Windows.Forms;
using OpenTK;
using Scene3D;

namespace _086shader
{
  public partial class Form1 : Form
  {
    static readonly string rev = "$Rev$".Split( ' ' )[ 1 ];

    /// <summary>
    /// Scene read from file.
    /// </summary>
    protected SceneBrep scene = new SceneBrep();

    /// <summary>
    /// Scene center point.
    /// </summary>
    protected Vector3 center = Vector3.Zero;

    /// <summary>
    /// Scene diameter.
    /// </summary>
    protected float diameter = 4.0f;

    /// <summary>
    /// GLControl guard flag.
    /// </summary>
    bool loaded = false;

    public Form1 ()
    {
      InitializeComponent();

      string param;
      Construction.InitParams( out param );
      textParam.Text = param ?? "";
      Text += " (rev: " + rev + ')';

      InitShaderRepository();
    }

    private void glControl1_Load ( object sender, EventArgs e )
    {
      InitOpenGL();
      SetupViewport();

      loaded = true;
      Application.Idle += new EventHandler( Application_Idle );
    }

    private void glControl1_Resize ( object sender, EventArgs e )
    {
      if ( !loaded ) return;

      SetupViewport();
      glControl1.Invalidate();
    }

    private void glControl1_Paint ( object sender, PaintEventArgs e )
    {
      Render();
    }

    private void checkVsync_CheckedChanged ( object sender, EventArgs e )
    {
      glControl1.VSync = checkVsync.Checked;
    }

    private void buttonOpen_Click ( object sender, EventArgs e )
    {
      OpenFileDialog ofd = new OpenFileDialog();

      ofd.Title = "Open Scene File";
      ofd.Filter = "Wavefront OBJ Files|*.obj;*.obj.gz" +
                   "|All scene types|*.obj";

      ofd.FilterIndex = 1;
      ofd.FileName = "";
      if ( ofd.ShowDialog() != DialogResult.OK )
        return;

      WavefrontObj objReader = new WavefrontObj();
      objReader.MirrorConversion = false;

      int faces = objReader.ReadBrep( ofd.FileName, scene );

      scene.BuildCornerTable();
      diameter = scene.GetDiameter( out center );
      scene.GenerateColors( 12 );
      ResetCamera();

      labelFile.Text = String.Format( "{0}: {1} faces", ofd.SafeFileName, faces );
      PrepareDataBuffers();
      glControl1.Invalidate();
    }

    private void buttonGenerate_Click ( object sender, EventArgs e )
    {
      Cursor.Current = Cursors.WaitCursor;

      scene.Reset();
      Construction cn = new Construction();

      int faces = cn.AddMesh( scene, Matrix4.Identity, textParam.Text );
      diameter = scene.GetDiameter( out center );

      if ( checkMulti.Checked )
      {
        Matrix4 translation, rotation;

        Matrix4.CreateTranslation( diameter, 0.0f, 0.0f, out translation );
        Matrix4.CreateRotationX( 90.0f, out rotation );
        faces += cn.AddMesh( scene, translation * rotation, textParam.Text );

        Matrix4.CreateTranslation( 0.0f, diameter, 0.0f, out translation );
        faces += cn.AddMesh( scene, translation, textParam.Text );

        Matrix4.CreateTranslation( diameter, diameter, 0.0f, out translation );
        faces += cn.AddMesh( scene, translation, textParam.Text );

        Matrix4.CreateTranslation( 0.0f, 0.0f, diameter, out translation );
        faces += cn.AddMesh( scene, translation, textParam.Text );

        Matrix4.CreateTranslation( diameter, 0.0f, diameter, out translation );
        faces += cn.AddMesh( scene, translation, textParam.Text );

        Matrix4.CreateTranslation( 0.0f, diameter, diameter, out translation );
        faces += cn.AddMesh( scene, translation, textParam.Text );

        Matrix4.CreateTranslation( diameter, diameter, diameter, out translation );
        faces += cn.AddMesh( scene, translation, textParam.Text );

        diameter = scene.GetDiameter( out center );
      }

      scene.BuildCornerTable();
      int errors = scene.CheckCornerTable( null );

      scene.GenerateColors( 12 );
      ResetCamera();

      labelFile.Text = String.Format( "{0} faces, {1} errors", faces, errors );
      PrepareDataBuffers();
      glControl1.Invalidate();

      Cursor.Current = Cursors.Default;
    }
  }
}