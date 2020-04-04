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


void transformvolume(double **Args)
{
    double *Isize_d, *mean, *A, *Iin, *Iout, *ThreadID;
    int Isize[3]={0,0,0};

    int x,y,z;

    // Location of pixel which will be come the current pixel
    double Tlocalx;
    double Tlocaly;
    double Tlocalz;
    
    // Linear interpolation variables
    int xBas[8], yBas[8], zBas[8];
    double perc[8];
    double xCom, yCom, zCom;
    double color[8]={0,0,0,0,0,0,0,0};
    
    // X,Y,Z coordinates of current pixel
    double xd,yd,zd;

    // Variables to store 1D index
    int indexI;
    
    // loop variables
    int p;
    
    // Multiple threads, one does the odd the other even indexes
    int offset;
    
    Isize_d=Args[0];
    mean=Args[1];
    A=Args[2];
    Iin=Args[3];
    Iout=Args[4];
    ThreadID=Args[5];
    
    Isize[0] = (int)Isize_d[0]; 
    Isize[1] = (int)Isize_d[1]; 
    Isize[2] = (int)Isize_d[2];
    
    if(ThreadID[0]==1) { offset = 0; }
    if(ThreadID[0]==2) { offset = 1; }
     
    // Loop through all image pixel coordinates
    for (z=offset; z<Isize[2]; z=z+2)
    {
        for (y=0; y<Isize[1]; y++)
        {
            for (x=0; x<Isize[0]; x++)
            {
                xd=x-mean[0]; yd=y-mean[1]; zd=z-mean[2];
                
                Tlocalx = mean[0] + A[0] * xd + A[1] *yd + A[2] *zd + A[3] * 1;
                Tlocaly = mean[1] + A[4] * xd + A[5] *yd + A[6] *zd + A[7] * 1;
                Tlocalz = mean[2] + A[8] * xd + A[9] *yd + A[10] *zd + A[11] * 1;
                                
                // Determine the coordinates of the pixel(s) which will be come the current pixel
                // (using linear interpolation)  
                xBas[0]=(int) floor(Tlocalx); yBas[0]=(int) floor(Tlocaly); zBas[0]=(int) floor(Tlocalz); 
                xBas[1]=xBas[0]+0;      yBas[1]=yBas[0]+0;      zBas[1]=zBas[0]+1;
                xBas[2]=xBas[0]+0;      yBas[2]=yBas[0]+1;      zBas[2]=zBas[0]+0;
                xBas[3]=xBas[0]+0;      yBas[3]=yBas[0]+1;      zBas[3]=zBas[0]+1;
                xBas[4]=xBas[0]+1;      yBas[4]=yBas[0]+0;      zBas[4]=zBas[0]+0;
                xBas[5]=xBas[0]+1;      yBas[5]=yBas[0]+0;      zBas[5]=zBas[0]+1;
                xBas[6]=xBas[0]+1;      yBas[6]=yBas[0]+1;      zBas[6]=zBas[0]+0;
                xBas[7]=xBas[0]+1;      yBas[7]=yBas[0]+1;      zBas[7]=zBas[0]+1;

                for (p=0; p<8; p++)
                {
                    if(xBas[p]>=0&&xBas[p]<Isize[0]&&yBas[p]>=0&&yBas[p]<Isize[1]&&zBas[p]>=0&&zBas[p]<Isize[2]) 
                    { 
                        indexI=mindex3(xBas[p],yBas[p],zBas[p],Isize[0],Isize[1]); 
                        color[p]=Iin[indexI]; 
                    } 
                    else { color[p]=0; }
                }

                // Linear interpolation constants (percentages)
                xCom=Tlocalx-floor(Tlocalx); yCom=Tlocaly-floor(Tlocaly);  zCom=Tlocalz-floor(Tlocalz);
                perc[0]=(1-xCom) * (1-yCom) * (1-zCom);
                perc[1]=(1-xCom) * (1-yCom) * zCom;
                perc[2]=(1-xCom) * yCom * (1-zCom);
                perc[3]=(1-xCom) * yCom * zCom;
                perc[4]=xCom * (1-yCom) * (1-zCom);
                perc[5]=xCom * (1-yCom) * zCom;
                perc[6]=xCom * yCom * (1-zCom);
                perc[7]=xCom * yCom * zCom;

                // Set the current pixel value
                indexI=mindex3(x,y,z,Isize[0],Isize[1]);
                Iout[indexI]=0;
                for (p=0; p<8; p++)
                {
                    Iout[indexI]+=color[p]*perc[p];
                }
            }
        }
    }   
    // Set the thread finished variables
    if(ThreadID[0]==1) { WaitForThread1 = 0; }
    if(ThreadID[0]==2) { WaitForThread2 = 0; }
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
    
    // Transformation matrix
    double A[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    // Size of input image
    double Isize_d[3]={0,0,0};
    const mwSize *dims;
    
    double mean[3]={0,0,0};
    
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
  Isize_d[0] = (double)dims[0]; Isize_d[1] = (double)dims[1]; Isize_d[2] = (double)dims[2];
  
  // Create output array
  plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);

  /* Assign pointers to each input. */
  Iin=mxGetPr(prhs[0]);
  M=mxGetPr(prhs[1]);
 
  A[0] = M[mindex2(0,0,4)];  A[1] = M[mindex2(0,1,4)];  A[2] = M[mindex2(0,2,4)];  A[3] = M[mindex2(0,3,4)];
  A[4] = M[mindex2(1,0,4)];  A[5] = M[mindex2(1,1,4)];  A[6] = M[mindex2(1,2,4)];  A[7] = M[mindex2(1,3,4)];
  A[8] = M[mindex2(2,0,4)];  A[9] = M[mindex2(2,1,4)];  A[10] = M[mindex2(2,2,4)]; A[11] = M[mindex2(2,3,4)];
  A[12] = M[mindex2(3,0,4)]; A[13] = M[mindex2(3,1,4)]; A[14] = M[mindex2(3,2,4)]; A[15] = M[mindex2(3,3,4)];
   
  
  /* Assign pointer to output. */
  Iout = mxGetPr(plhs[0]);
  
  /* Center of the volume */
  mean[0]=Isize_d[0]/2;  mean[1]=Isize_d[1]/2;  mean[2]=Isize_d[2]/2;
  
  WaitForThread1 = 1;
  ThreadArgs1[0]=Isize_d;
  ThreadArgs1[1]=mean;
  ThreadArgs1[2]=A;
  ThreadArgs1[3]=Iin;
  ThreadArgs1[4]=Iout;
  ThreadArgs1[5]=ThreadID1;
  _beginthread( transformvolume, 0, ThreadArgs1 );
  
  WaitForThread2 = 1;
  ThreadArgs2[0]=Isize_d;
  ThreadArgs2[1]=mean;
  ThreadArgs2[2]=A;
  ThreadArgs2[3]=Iin;
  ThreadArgs2[4]=Iout;
  ThreadArgs2[5]=ThreadID2;
  _beginthread( transformvolume, 0, ThreadArgs2 );
  
  // Wait for the two threads to finish
  while( WaitForThread1 )  {  Sleep( 5 ); }
  while( WaitForThread2 )  {  Sleep( 5 ); }
}
        

