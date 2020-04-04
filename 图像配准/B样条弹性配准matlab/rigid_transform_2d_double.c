#include "mex.h"
#include "math.h"
// undef needed for LCC compiler
#undef EXTERN_C
#include <windows.h>
#include <process.h>    

/* Imrigide transform Iout=rigide_transform_3d_double(Iin,Minv)
  * This function transforms a volume with a 4x4 transformation matrix 
  * 
  * Function is written by D.Kroon University of Twente (September 2008)
  */

// Variables used to detect if the threads are finished
// (volatile: reload the variable instead of using the value available in a register)
static volatile int WaitForThread1;
static volatile int WaitForThread2;

// Convert 2D/3D matrix index to 1D index
int mindex2(int x, int y, int sizx) { return y*sizx+x; }
int mindex3(int x, int y, int z, int sizx, int sizy) { return z*sizx*sizy+y*sizx+x;}


unsigned __stdcall  transformvolume(double **Args)
{
    double *Isize_d, *mean, *A, *Iin, *Iout, *ThreadID;
    int Isize[2]={0,0};

    int x,y;

    // Location of pixel which will be come the current pixel
    double Tlocalx;
    double Tlocaly;
    
    // Linear interpolation variables
    int xBas[4], yBas[4];
    double perc[4];
    double xCom, yCom;
    double color[4]={0,0,0,0};
    
    // X,Y,Z coordinates of current pixel
    double xd,yd;

    // Variables to store 1D index
    int indexI;
    
    // Multiple threads, one does the odd the other even indexes
    int offset;
    //int start;
    //int end;
    
    Isize_d=Args[0];
    mean=Args[1];
    A=Args[2];
    Iin=Args[3];
    Iout=Args[4];
    ThreadID=Args[5];
    
    Isize[0] = (int)Isize_d[0]; 
    Isize[1] = (int)Isize_d[1]; 
    
    if(ThreadID[0]==1) { offset=0; } //start=0; end = Isize[1]/2; 
    if(ThreadID[0]==2) { offset=1; } //start=Isize[1]/2; end=Isize[1];
     
    // Loop through all image pixel coordinates
    for (y=offset; y<Isize[1]; y=y+2)
    {
        for (x=0; x<Isize[0]; x++)
        {
            xd=(double)x-mean[0]; yd=(double)y-mean[1];

            Tlocalx = mean[0] + A[0] * xd + A[1] *yd + A[2] * 1;
            Tlocaly = mean[1] + A[3] * xd + A[4] *yd + A[5] * 1;

            // Determine the coordinates of the pixel(s) which will be come the current pixel
            // (using linear interpolation)  
            xBas[0]=(int) floor(Tlocalx); yBas[0]=(int) floor(Tlocaly);
            xBas[1]=xBas[0]+0;      yBas[1]=yBas[0]+1;
            xBas[2]=xBas[0]+1;      yBas[2]=yBas[0]+0;
            xBas[3]=xBas[0]+1;      yBas[3]=yBas[0]+1;

            if(xBas[0]>=0&&xBas[0]<Isize[0]&&yBas[0]>=0&&yBas[0]<Isize[1]) { indexI=mindex2(xBas[0],yBas[0],Isize[0]); color[0]=Iin[indexI]; } else { color[0]=0;}
            if(xBas[1]>=0&&xBas[1]<Isize[0]&&yBas[1]>=0&&yBas[1]<Isize[1]) { indexI=mindex2(xBas[1],yBas[1],Isize[0]); color[1]=Iin[indexI]; } else { color[1]=0;}
            if(xBas[2]>=0&&xBas[2]<Isize[0]&&yBas[2]>=0&&yBas[2]<Isize[1]) { indexI=mindex2(xBas[2],yBas[2],Isize[0]); color[2]=Iin[indexI]; } else { color[2]=0;}
            if(xBas[3]>=0&&xBas[3]<Isize[0]&&yBas[3]>=0&&yBas[3]<Isize[1]) { indexI=mindex2(xBas[3],yBas[3],Isize[0]); color[3]=Iin[indexI]; } else { color[3]=0;}

            // Linear interpolation constants (percentages)
            xCom=Tlocalx-floor(Tlocalx); yCom=Tlocaly-floor(Tlocaly);
            perc[0]=(1-xCom) * (1-yCom);
            perc[1]=(1-xCom) * yCom;
            perc[2]=xCom * (1-yCom);
            perc[3]=xCom * yCom;

            // Set the current pixel value
            indexI=mindex2(x,y,Isize[0]);
            Iout[indexI]=color[0]*perc[0]+color[1]*perc[1]+color[2]*perc[2]+color[3]*perc[3];
        }
    }
    // Set the thread finished variables
    if(ThreadID[0]==1) { WaitForThread1 = 0; }
    if(ThreadID[0]==2) { WaitForThread2 = 0; }
    
    // explicit end thread, helps to ensure proper recovery of resources allocated for the thread
    _endthreadex( 0 );
    return 0;
    
}

// The matlab mex function
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    // Ox and Oy are the grid points
    // Zo is the input image
    // Zi is the transformed image
    // nx and ny are the number of grid points (inside the image)
    double *Iin, *Iout, *M;
    
    // double pointer array to store all needed function variables
    double **ThreadArgs1,**ThreadArgs2;
    
    double ThreadID1[1]={1}; // ID of first Thread 
    double ThreadID2[1]={2}; // ID of second Thread
    
    HANDLE *ThreadList; // Handles to the worker threads
    
    // Transformation matrix
    double A[9]={0,0,0,0,0,0,0,0,0};

    // Size of input image
    double Isize_d[2]={0,0};
    const mwSize *dims;
    
    double mean[2]={0,0};
    
    // Reserve room for handles of threads in ThreadList
    ThreadList = (HANDLE*)malloc(2* sizeof( HANDLE ));
    
     // Reserve room for 6 function variables(arrays)
    ThreadArgs1 = (double **)malloc( 6* sizeof( double * ) );  
    ThreadArgs2 = (double **)malloc( 6* sizeof( double * ) );  
    
    
  /* Check for proper number of arguments. */
  if(nrhs!=2) {
    mexErrMsgTxt("Two inputs are required.");
  } else if(nlhs!=1) {
    mexErrMsgTxt("One output required");
  }
  // nsubs=mxGetNumberOfDimensions(prhs[0]);
     
  // Get the sizes of the image
  dims = mxGetDimensions(prhs[0]);   
  Isize_d[0] = (double)dims[0]; Isize_d[1] = (double)dims[1]; 
  
  // Create output array
  plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);

  /* Assign pointers to each input. */
  Iin=mxGetPr(prhs[0]);
  M=mxGetPr(prhs[1]);
 
  A[0] = M[mindex2(0,0,3)]; A[1] = M[mindex2(0,1,3)]; A[2] = M[mindex2(0,2,3)]; 
  A[3] = M[mindex2(1,0,3)]; A[4] = M[mindex2(1,1,3)]; A[5] = M[mindex2(1,2,3)]; 
  A[6] = M[mindex2(2,0,3)]; A[7] = M[mindex2(2,1,3)]; A[8] = M[mindex2(2,2,3)]; 
  
  
  /* Assign pointer to output. */
  Iout = mxGetPr(plhs[0]);
  
  /* Center of the volume */
  mean[0]=Isize_d[0]/2;  mean[1]=Isize_d[1]/2;  
  
  WaitForThread1 = 1;
  ThreadArgs1[0]=Isize_d;
  ThreadArgs1[1]=mean;
  ThreadArgs1[2]=A;
  ThreadArgs1[3]=Iin;
  ThreadArgs1[4]=Iout;
  ThreadArgs1[5]=ThreadID1;
  ThreadList[0] = (HANDLE)_beginthreadex( NULL, 0, &transformvolume, ThreadArgs1 , 0, NULL );
  
  WaitForThread2 = 1;
  ThreadArgs2[0]=Isize_d;
  ThreadArgs2[1]=mean;
  ThreadArgs2[2]=A;
  ThreadArgs2[3]=Iin;
  ThreadArgs2[4]=Iout;
  ThreadArgs2[5]=ThreadID2;
  ThreadList[1] = (HANDLE)_beginthreadex( NULL, 0, &transformvolume, ThreadArgs2 , 0, NULL );
  
  WaitForSingleObject(ThreadList[0], INFINITE);
  WaitForSingleObject(ThreadList[1], INFINITE);
  
  CloseHandle( ThreadList[0] );
  CloseHandle( ThreadList[1] );
}
        

