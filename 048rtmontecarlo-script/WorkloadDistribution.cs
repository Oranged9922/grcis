﻿using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Net;
using System.Net.Sockets;
using System.Text;
using System.Threading;
using MathSupport;
using _048rtmontecarlo;

namespace Rendering
{
  /// <summary>
  /// Takes care of distribution of rendering work between local threads and remote/network render clients
  /// </summary>
  public class Master
  {
    public static Master instance; // singleton

    public ConcurrentQueue<Assignment> availableAssignments;

    private List<NetworkWorker> networkWorkers;

    private Thread[] pool;

    public Thread mainThread;

    public int totalNumberOfAssignments;

    public int finishedAssignments;

    // width and height of one block of pixels (rendered at one thread at the time); 64 seems to be optimal; should be power of 2 and larger than 8
    public const int assignmentSize = 64;

    public Progress progressData;

    public int assignmentRoundsFinished = 0;
    public int assignmentRoundsTotal;

    public Bitmap    bitmap;
    public IRayScene scene;
    public IRenderer renderer;

    /// <summary>
    /// Constructor which takes also care of initializing assignments
    /// </summary>
    /// <param name="bitmap">Main bitmap - used also in PictureBox in Form1</param>
    /// <param name="scene">Scene to render</param>
    /// <param name="renderer">Rendered to use for RenderPixel method</param>
    public Master ( Bitmap bitmap, IRayScene scene, IRenderer renderer )
    {
      finishedAssignments = 0;
   
      this.scene    = scene;
      this.renderer = renderer;
      this.bitmap = bitmap;     

      if ( RenderClientsForm.instance == null )
      {
        RenderClientsForm.instance = new RenderClientsForm ();
      }
    }

    /// <summary>
    /// Creates threadpool and starts all threads on Consume method
    /// </summary>
    /// <param name="threads">Number of threads to be used for rendering</param>
    public void StartThreads ( int threads )
    {
      pool = new Thread[threads];

      AssignNetworkWorkerToStream ();

      Thread imageFetcher = new Thread ( RenderedImageReceiver );
      imageFetcher.Name = "ImageFetcher";
      imageFetcher.Start ();

      for ( int i = 0; i < threads; i++ )
      {
        Thread newThread = new Thread ( Consume );
        newThread.Name = "RenderThread " + i;
        newThread.Priority = ThreadPriority.AboveNormal;
        pool [ i ]         = newThread;
        newThread.Start ();
      }

      mainThread = pool [ 0 ];

      for ( int i = 0; i < (int) threads; i++ )
      {
        pool [ i ].Join ();
        pool [ i ] = null;
      }

      imageFetcher.Join ();
      imageFetcher = null;
    }

    /// <summary>
    /// Consumer-producer based multithreading work distribution
    /// Each thread waits for a new Assignment to be added to availableAssignments queue
    /// Most of the time is number of items in availableAssignments expected to be several times larger than number of threads
    /// </summary>
    protected void Consume ()
    {
      MT.InitThreadData ();

      while ( finishedAssignments < totalNumberOfAssignments )
      {
        Assignment newAssignment;
        availableAssignments.TryDequeue ( out newAssignment );

        if ( !progressData.Continue ) // test whether rendering should end (Stop button pressed) 
          return;

        if ( newAssignment == null ) // TryDequeue was not successful
          continue;

        float[] colorArray = newAssignment.Render ( false, renderer, progressData );
        BitmapMerger ( colorArray, newAssignment.x1, newAssignment.y1, newAssignment.x2 + 1, newAssignment.y2 + 1 );

        if ( newAssignment.stride == 1 )
        {
          finishedAssignments++;
          assignmentRoundsFinished++;
        }
        else
        {
          newAssignment.stride = newAssignment.stride >> 1; // stride values: 8 > 4 > 2 > 1
          assignmentRoundsFinished++;
          availableAssignments.Enqueue ( newAssignment );
        }
      }
    }

    /// <summary>
    /// Creates new assignments based on width and heigh of bitmap and assignmentSize
    /// </summary>
    /// <param name="bitmap">Main bitmap - used also in PictureBox in Form1</param>
    /// <param name="scene">Scene to render</param>
    /// <param name="renderer">Rendered to use for RenderPixel method</param>
    public void InitializeAssignments ( Bitmap bitmap, IRayScene scene, IRenderer renderer )
    {
      availableAssignments = new ConcurrentQueue<Assignment> ();

      int numberOfAssignmentsOnWidth = bitmap.Width % assignmentSize == 0
        ? bitmap.Width / assignmentSize
        : bitmap.Width / assignmentSize + 1;

      int numberOfAssignmentsOnHeight = bitmap.Height % assignmentSize == 0
        ? bitmap.Height / assignmentSize
        : bitmap.Height / assignmentSize + 1;


      for ( int y = 0; y < numberOfAssignmentsOnHeight; y++ )
      {
        for ( int x = 0; x < numberOfAssignmentsOnWidth; x++ )
        {
          int localX = x * assignmentSize;
          int localY = y * assignmentSize;

          Assignment newAssignment = new Assignment ( localX,
                                                      localY,
                                                      localX + assignmentSize - 1,
                                                      localY + assignmentSize - 1,
                                                      bitmap.Width, 
                                                      bitmap.Height );
          availableAssignments.Enqueue ( newAssignment );
        }
      }


      totalNumberOfAssignments = availableAssignments.Count;
      assignmentRoundsTotal    = totalNumberOfAssignments * 4;
    }

    /// <summary>
    /// Goes through all clients from RenderClientsForm and assigns a NetworkWorker to each of them
    /// </summary>
    public void AssignNetworkWorkerToStream ()
    {
      foreach ( Client client in RenderClientsForm.instance.clients )
      {
        NetworkWorker newWorker = new NetworkWorker ( client.address );
        networkWorkers = new List<NetworkWorker> ();

        if ( !newWorker.ConnectToClient () ) // removes NetworkWorker instance in case of failure to connect to the client
        {
          newWorker = null;
        }
        else
        {
          networkWorkers.Add ( newWorker );
          newWorker.SendNecessaryObjects ();
          newWorker.TryToGetNewAssignment ();
          newWorker.TryToGetNewAssignment ();
          newWorker.TryToGetNewAssignment ();
          newWorker.TryToGetNewAssignment ();
          newWorker.TryToGetNewAssignment ();

          newWorker.SendEndingAssignment ();
        }
      }
    }

    /// <summary>
    /// Adds colors represented in newBitmap array to main bitmap
    /// </summary>
    /// <param name="newBitmap">Float values (for later HDR support) representing pixel color values</param>
    public void BitmapMerger ( float[] newBitmap, int x1, int y1, int x2, int y2 )
    {
      lock ( bitmap )
      {
        int arrayPosition = 0;

        for ( int y = y1; y < Math.Min ( y2, bitmap.Height ); y++ )
        {
          for ( int x = x1; x < Math.Min ( x2, bitmap.Width ); x++ )
          {
            if ( !float.IsInfinity ( newBitmap[arrayPosition] ) )
            {
              Color color = Color.FromArgb ( Math.Min ( (int) newBitmap [ arrayPosition ], 255 ),
                                             Math.Min ( (int) newBitmap [ arrayPosition + 1 ], 255 ),
                                             Math.Min ( (int) newBitmap [ arrayPosition + 2 ], 255 ) );
              bitmap.SetPixel ( x, y, color );
            }

            arrayPosition += 3;
          }
        }
      }
    }

    /// <summary>
    /// Goes through NetworkStreams of all clients and checks if there is rendered image pending in them
    /// Plus it takes care of that data via ReceiveRenderedImage method
    /// </summary>
    public void RenderedImageReceiver ()
    {
      if ( networkWorkers == null || networkWorkers.Count == 0 )
      {
        return;
      }

      while ( finishedAssignments < totalNumberOfAssignments )
      {
        foreach ( NetworkWorker worker in networkWorkers )
        {
          worker.ReceiveRenderedImage ();
        }
      }
    }
  }


  /// <summary>
  /// Takes care of network communication with with 1 render client
  /// </summary>
  class NetworkWorker
  {
    private readonly IPAddress  ipAdr;
    private          IPEndPoint endPoint;
    private const    int        port = 5000;

    private TcpClient     client;
    private NetworkStream stream;

    public NetworkWorker ( IPAddress ipAdr )
    {
      this.ipAdr = ipAdr;
    }

    /// <summary>
    /// Establishes NetworkStream with desired client
    /// </summary>
    /// <returns>True if connection was succesfull, False otherwise</returns>
    public bool ConnectToClient ()
    {
      if ( ipAdr == null )
        return false;

      if ( endPoint == null )
        endPoint = new IPEndPoint ( ipAdr, port );

      client = new TcpClient ();

      try
      {
        client.Connect ( endPoint );
      }
      catch ( SocketException e )
      {
        return false;
      }

      stream = client.GetStream ();

      client.ReceiveBufferSize = 1024 * 1024 * 32;
      client.SendBufferSize    = 1024 * 1024 * 32;
      //client.NoDelay = true;

      return true;
    }


    /// <summary>
    /// Sends all objects which are necessary to render scene to client
    /// IRayScene - scene representation like solids, textures, lights, camera, ...
    /// IRenderer - renderer itself including IImageFunction; needed for RenderPixel method
    /// </summary>
    public void SendNecessaryObjects ()
    {
      //lock ( stream )
      {
        NetworkSupport.SendObject<IRayScene> ( Master.instance.scene, client, stream );
      }

      //lock ( stream )
      {
        NetworkSupport.SendObject<IRenderer> ( Master.instance.renderer, client, stream );
      }      
    }


    private const int bufferSize = ( Master.assignmentSize * Master.assignmentSize * 3 + 2) * sizeof ( float );
    /// <summary>
    /// Checks whether there is any data in NetworkStream to read
    /// If so, it reads it as - expected format is array of floats
    ///   - first 2 floats represents x1 and y1 coordinates - position in main bitmap;
    ///   - rest of array are colors of rendered bitmap - 3 floats (RGB values) per pixel;
    ///   - stored per lines from left upper corner (coordinates position)
    /// </summary>
    public void ReceiveRenderedImage ()
    {
      if ( stream.DataAvailable )
      {
        int totalReceivedSize = 0;
        int leftToReceive = bufferSize;

        byte[] receiveBuffer = new byte[bufferSize];        

        while ( leftToReceive > 0 )  // Loop until enough data is received
        {
          int latestReceivedSize = stream.Read ( receiveBuffer, totalReceivedSize, leftToReceive );
          leftToReceive -= latestReceivedSize;
          totalReceivedSize += latestReceivedSize;
        }

        // Use parts of bigBuffer - separate and convert data to coordinates and floats representing colors of pixels
        float[] coordinates = new float[2];
        float[] floatBuffer = new float[Master.assignmentSize * Master.assignmentSize * 3];
        Buffer.BlockCopy ( receiveBuffer, 0, coordinates, 0, 2 * sizeof ( float ) );
        Buffer.BlockCopy ( receiveBuffer, 2 * sizeof ( float ), floatBuffer, 0, floatBuffer.Length * sizeof ( float ) );

        Master.instance.BitmapMerger ( floatBuffer, 
                                       (int) coordinates [ 0 ],
                                       (int) coordinates [ 1 ],
                                       (int) coordinates [ 0 ] + Master.assignmentSize,
                                       (int) coordinates [ 1 ] + Master.assignmentSize );

        Master.instance.finishedAssignments++;
        //TryToGetNewAssignment();
      }
    }    

    public void TryToGetNewAssignment ()
    {
      Assignment newAssignment = null;

      while ( Master.instance.finishedAssignments < Master.instance.totalNumberOfAssignments - 7 )
      {
        Master.instance.availableAssignments.TryDequeue ( out newAssignment );

        if ( !Master.instance.progressData.Continue ) // test whether rendering should end (Stop button pressed) 
          return;

        if ( newAssignment == null ) // TryDequeue was not succesfull
          continue;

        lock ( stream )
        {
          NetworkSupport.SendObject<Assignment> ( newAssignment, client, stream );
        }
        
        break;
      }
    }

    /// <summary>
    /// Send dummy Assignment to signal client that rendering ended and/or client should not wait for more assignments
    /// </summary>
    public void SendEndingAssignment ()
    {
      Assignment newAssignment = new Assignment ( -1, -1, -1, -1, -1, -1 );

      lock ( stream )
      {
        NetworkSupport.SendObject<Assignment> ( newAssignment, client, stream );
      }
    }
  }


  /// <summary>
  /// Represents 1 render work (= square of pixels to render at specific stride)
  /// </summary>
  [Serializable]
  public class Assignment
  {
    public int x1, y1, x2, y2;

    public const int assignmentSize = 64;

    public int stride; // stride of 'n' means that only each 'n'-th pixel is rendered (for sake of dynamic rendering)

    private readonly int bitmapWidth, bitmapHeight;

    public Assignment ( int x1, int y1, int x2, int y2, int bitmapWidth, int bitmapHeight )
    {
      this.x1 = x1;
      this.y1 = y1;
      this.x2 = x2;
      this.y2 = y2;
      this.bitmapWidth = bitmapWidth;
      this.bitmapHeight = bitmapHeight;

      // stride values: 8 > 4 > 2 > 1; initially always 8
      // decreases at the end of rendering of current assignment and therefore makes another render of this assignment more detailed
      stride = 8;   
    }

    /// <summary>
    /// Main render method
    /// Directly writes pixel colors to the main bitmap after rendering them
    /// </summary>
    public float[] Render ( bool renderEverything, IRenderer renderer, Progress progressData = null )
    {
      float[] returnArray = new float[assignmentSize * assignmentSize * 3];

      if ( renderEverything )
      {
        stride = 1;
      }

      for ( int y = y1; y <= y2; y += stride )
      {
        for ( int x = x1; x <= x2; x += stride )
        {
          if ( stride != 8 && ( y % ( stride << 1 ) == 0 ) && ( x % ( stride << 1 ) == 0 ) && !renderEverything) // prevents rendering of already rendered pixels
          {
            returnArray[PositionInArray ( x, y )]     = float.PositiveInfinity;
            returnArray[PositionInArray ( x, y ) + 1] = float.PositiveInfinity;
            returnArray[PositionInArray ( x, y ) + 2] = float.PositiveInfinity;
            continue;
          }

          if ( x >= bitmapWidth || y >= bitmapHeight )
          {
            continue;
          }

          double[] color = new double[3];

          renderer.RenderPixel ( x, y, color ); // called at desired IRenderer; gets pixel color

          if ( stride == 1 )
          {
            returnArray [ PositionInArray ( x, y ) ] = (float) ( color [ 0 ] * 255.0 );
            returnArray [ PositionInArray ( x, y ) + 1 ] = (float) ( color [ 1 ] * 255.0 );
            returnArray [ PositionInArray ( x, y ) + 2 ] = (float) ( color [ 2 ] * 255.0 );
          }
          else
          {
            for ( int j = y; j < y + stride; j++ )
            {
              if ( j < bitmapHeight )
              {
                for ( int i = x; i < x + stride; i++ )
                {
                  if ( i < bitmapWidth )
                  {
                    // actual set of pixel color to main bitmap
                    returnArray [ PositionInArray ( i, j ) ] = (float) ( color [ 0 ] * 255.0 );
                    returnArray [ PositionInArray ( i, j ) + 1 ] = (float) ( color [ 1 ] * 255.0 );
                    returnArray [ PositionInArray ( i, j ) + 2 ] = (float) ( color [ 2 ] * 255.0 );
                  }
                }
              }
            }
          }

          if ( progressData != null )
          {
            lock ( Master.instance.progressData )
            {
              // test whether rendering should end (Stop button pressed) 
              if ( !Master.instance.progressData.Continue )
                return returnArray;

              // synchronization of bitmap with PictureBox in Form and update of progress (percentage of done work)
              if ( Master.instance.mainThread == Thread.CurrentThread )
              {
                Master.instance.progressData.Finished =
                  Master.instance.assignmentRoundsFinished / (float) Master.instance.assignmentRoundsTotal;
                Master.instance.progressData.Sync ( Master.instance.bitmap );
              }
            }
          }         
        }
      }

      return returnArray;
    }

    /// <summary>
    /// Computes position of coordinates (bitmap) into one-dimensional array of floats (with 3 floats (RGB channels) per pixel)
    /// </summary>
    /// <param name="x">X coordinate in bitmap</param>
    /// <param name="y">Y coordinate in bitmap</param>
    /// <returns>Position in array of floats</returns>
    private int PositionInArray ( int x, int y )
    {
      return ( ( y - y1 ) * assignmentSize + ( x - x1 ) ) * 3;
    }
  }
}