#include "mex.h"
#include "math.h"
// undef needed for LCC compiler
#undef EXTERN_C
#include <windows.h>
#include <process.h>    

/* 3D Bspline transformation grid function
 * function [Vout,Tx,Ty,Tz]t=bspline_transform_3d_double(Ox,Oy,Oz,Vin,dx,dy,dz)
 * 
 * Ox, Oy, Oz are the grid points coordinates
 * Vin is input image, Vout the transformed output image
 * dx, dy and dz are the spacing of the b-spline knots
 *
 * Iout: The transformed image
 * Tx: The transformation field in x direction
 * Ty: The transformation field in y direction
 * Tz: The transformation field in y direction
 *
 * This function is an implementation of the b-spline registration
 * algorithm in "D. Rueckert et al. : Nonrigid Registration Using Free-Form 
 * Deformations: Application to Breast MR Images".
 * 
 * We used "Fumihiko Ino et al. : a data distrubted parallel algortihm for 
 * nonrigid image registration" for the correct formula's, because 
 * (most) other papers contain errors. 
 *
 * Function is written by D.Kroon University of Twente (July 2008)
 */

// Variables used to detect if the threads are finished
// (volatile: reload the variable instead of using the value available in a register)
static volatile int WaitForThread1;
static volatile int WaitForThread2;

// Convert 2D/3D matrix index to 1D index
int mindex2(int x, int y, int sizx) { return y*sizx+x; }
int mindex3(int x, int y, int z, int sizx, int sizy) { return z*sizx*sizy+y*sizx+x;}


unsigned __stdcall transformvolume(double **Args)
{
    double *Bu, *Bv, *Bw, *ZI, *Tx, *Ty, *Tz;
    double *dxa, *dya, *dza, *ThreadID, *Ox, *Oy, *Oz, *ZO;
    double *ZOsize_d;
    double *Osize_d;
    double *nlhs_d;
    int ZOsize[3]={0,0,0};
    int Osize[3]={0,0,0};
    
    // Multiple threads, one does the odd the other even indexes
    int offset;
        
    // Location of pixel which will be come the current pixel
    double Tlocalx;
    double Tlocaly;
    double Tlocalz;
    
    
    // Variables to store 1D index
    int indexO;
    int indexZI;
    int indexZO;
    
    // Grid distance;
    int dx,dy,dz; 
    
    // X,Y,Z coordinates of current pixel
    int x,y,z;
    
    // Linear interpolation variables
    // Linear interpolation variables
    int xBas[8], yBas[8], zBas[8];
    double perc[8];
    double xCom, yCom, zCom;
    double color[8]={0,0,0,0,0,0,0,0};
    
    // B-spline variables
    int u_index=0; 
    int v_index=0;
    int w_index=0;
    
    int i, j, k;
    
    // B-Spline loop variabels
    int l,m,n,p;
    int nlhs=0;
    
    // Split input into variables
    Bu=Args[0];
    Bv=Args[1];
    Bw=Args[2];
    ZOsize_d=Args[3];
    Osize_d=Args[4];
    ZI=Args[5];
    Tx=Args[6];
    Ty=Args[7];
    Tz=Args[8];
    dxa=Args[9];
    dya=Args[10];
    dza=Args[11];
    ThreadID=Args[12];
    Ox=Args[13];
    Oy=Args[14];
    Oz=Args[15];
    ZO=Args[16];
    nlhs_d=Args[17];
   
    nlhs=(int)nlhs_d[0];
    ZOsize[0] = (int)ZOsize_d[0]; 
    ZOsize[1] = (int)ZOsize_d[1]; 
    ZOsize[2] = (int)ZOsize_d[2]; 
    
    Osize[0] = (int)Osize_d[0]; 
    Osize[1] = (int)Osize_d[1]; 
    Osize[2] = (int)Osize_d[1]; 
    
    /* Get the spacing of the uniform b-spline grid */
    dx=(int)dxa[0]; dy=(int)dya[0]; dz=(int)dza[0];
  
    
    if(ThreadID[0]==1) { offset = 0; }
    if(ThreadID[0]==2) { offset = 1; }
    
    
 // Loop through all image pixel coordinates
     for (z=offset; z<ZOsize[2]; z=z+2)
     {
        for (y=0; y<ZOsize[1]; y++)
        {
            for (x=0; x<ZOsize[0]; x++)
            {
                // Calculate the indexes need to loop up the B-spline values.
                u_index=x%dx; 
                v_index=y%dy;
                w_index=z%dz;

                i=(int)floor(x/dx); // (first row outside image against boundary artefacts)
                j=(int)floor(y/dy);
                k=(int)floor(z/dz); 
 
                // This part calculates the coordinates of the pixel
                // which will be transformed to the current x,y pixel.
                Tlocalx=0; Tlocaly=0; Tlocalz=0;
                for(l=0; l<4; l++)
                {
                    for(m=0; m<4; m++)
                    {    
                        for(n=0; n<4; n++)
                        {       
                            if(((i+l)>=0)&&((i+l)<Osize[0])&&((j+m)>=0)&&((j+m)<Osize[1])&&((k+n)>=0)&&((k+n)<Osize[2]))
                            {
                                 indexO=mindex3(i+l,j+m,k+n,Osize[0],Osize[1]);
                                 Tlocalx=Tlocalx+Bu[mindex2(l,u_index,4)]*Bv[mindex2(m,v_index,4)]*Bw[mindex2(n,w_index,4)]*Ox[indexO];
                                 Tlocaly=Tlocaly+Bu[mindex2(l,u_index,4)]*Bv[mindex2(m,v_index,4)]*Bw[mindex2(n,w_index,4)]*Oy[indexO];
                                 Tlocalz=Tlocalz+Bu[mindex2(l,u_index,4)]*Bv[mindex2(m,v_index,4)]*Bw[mindex2(n,w_index,4)]*Oz[indexO];
                            }
                        }
                    }
                }            

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
                    if(xBas[p]>=0&&xBas[p]<ZOsize[0]&&yBas[p]>=0&&yBas[p]<ZOsize[1]&&zBas[p]>=0&&zBas[p]<ZOsize[2]) 
                    { 
                        indexZI=mindex3(xBas[p],yBas[p],zBas[p],ZOsize[0],ZOsize[1]); 
                        color[p]=ZO[indexZI]; 
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
                indexZO=mindex3(x,y,z,ZOsize[0],ZOsize[1]);
                ZI[indexZO]=0;
                for (p=0; p<8; p++)
                {
                    ZI[indexZO]+=color[p]*perc[p];
                }
               
                // Store transformation field
                if(nlhs>1) { Tx[indexZO]=Tlocalx-(double)x; }
                if(nlhs>2) { Ty[indexZO]=Tlocaly-(double)y; }
                if(nlhs>3) { Tz[indexZO]=Tlocalz-(double)z; }
                
            }
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
    double *Ox,*Oy,*Oz,*ZO, *dxa, *dya,*dza,*ZI, *Tx,*Ty,*Tz;
    
        
    // double pointer array to store all needed function variables
    double **ThreadArgs1,**ThreadArgs2;
    
    HANDLE *ThreadList; // Handles to the worker threads
    
    double ThreadID1[1]={1}; // ID of first Thread 
    double ThreadID2[1]={2}; // ID of second Thread
    
    double nlhs_d[1]={0};
    
    // Size of input image
    mwSize  ZOsizex, ZOsizey, ZOsizez;
    double ZOsize_d[3]={0,0,0};
    const mwSize *dims;
    // Size of grid
    mwSize  Osizex, Osizey, Osizez;
    double Osize_d[3]={0,0,0};
   
    // B-spline variables
    double u,v,w;
    int u_index=0; 
    int v_index=0;
    int w_index=0;
    
    double *Bu, *Bv, *Bw;
    
    // Grid distance;
    int dx,dy,dz; 
    
    // X,Y,Z coordinates of current pixel
    int x,y,z;
    
    
  // Reserve room for 6 function variables(arrays)
  ThreadArgs1 = (double **)malloc( 18* sizeof( double * ) );  
  ThreadArgs2 = (double **)malloc( 18* sizeof( double * ) );  
    
  // Reserve room for handles of threads in ThreadList
  ThreadList = (HANDLE*)malloc(2* sizeof( HANDLE ));
    
  
  /* Check for proper number of arguments. */
  if(nrhs!=7) {
    mexErrMsgTxt("Seven inputs are required.");
  }
 
  // Get the sizes of the grid
  dims = mxGetDimensions(prhs[0]);   
  Osizex = dims[0]; 
  Osizey = dims[1];
  Osizez = dims[2];
  
  // Create image matrix for the return arguments with the size of input image   
  dims = mxGetDimensions(prhs[3]);  
  ZOsizex = dims[0]; 
  ZOsizey = dims[1];
  ZOsizez = dims[2];
  
  plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); 
  if(nlhs>1) { plhs[1] = mxCreateDoubleMatrix(ZOsizex,ZOsizey, mxREAL); }
  if(nlhs>2) { plhs[2] = mxCreateDoubleMatrix(ZOsizex,ZOsizey, mxREAL); }
  if(nlhs>3) { plhs[3] = mxCreateDoubleMatrix(ZOsizex,ZOsizey, mxREAL); }
  
  
  /* Assign pointers to each input. */
  Ox=mxGetPr(prhs[0]);
  Oy=mxGetPr(prhs[1]);
  Oz=mxGetPr(prhs[2]);
  ZO=mxGetPr(prhs[3]);
  dxa=mxGetPr(prhs[4]);
  dya=mxGetPr(prhs[5]);
  dza=mxGetPr(prhs[6]);
    
  /* Get the spacing of the uniform b-spline grid */
  dx=(int)dxa[0]; dy=(int)dya[0]; dz=(int)dza[0]; 
  

  /* Assign pointer to output. */
  ZI = mxGetPr(plhs[0]);
  if(nlhs>1) { Tx = mxGetPr(plhs[1]); }
  if(nlhs>2) { Ty = mxGetPr(plhs[2]); }
  if(nlhs>3) { Tz = mxGetPr(plhs[3]); }
  
   // Make polynomial look up tables 
  Bu=malloc(dx*4*sizeof(double));
  Bv=malloc(dy*4*sizeof(double));
  Bw=malloc(dz*4*sizeof(double));
  for (x=0; x<dx; x++)
  {
    u=((double)x/(double)dx)-floor((double)x/(double)dx);
    Bu[mindex2(0,x,4)] = pow((1-u),3)/6;
    Bu[mindex2(1,x,4)] = ( 3*pow(u,3) - 6*pow(u,2) + 4)/6;
    Bu[mindex2(2,x,4)] = (-3*pow(u,3) + 3*pow(u,2) + 3*u + 1)/6;
    Bu[mindex2(3,x,4)] = pow(u,3)/6;
  }
  
  for (y=0; y<dy; y++)
  {
    v=((double)y/(double)dy)-floor((double)y/(double)dy);
    Bv[mindex2(0,y,4)] = pow((1-v),3)/6;
    Bv[mindex2(1,y,4)] = ( 3*pow(v,3) - 6*pow(v,2) + 4)/6;
    Bv[mindex2(2,y,4)] = (-3*pow(v,3) + 3*pow(v,2) + 3*v + 1)/6;
    Bv[mindex2(3,y,4)] = pow(v,3)/6;
  }
  
  for (z=0; z<dz; z++)
  {
    w=((double)z/(double)dz)-floor((double)z/(double)dz);
    Bw[mindex2(0,z,4)] = pow((1-w),3)/6;
    Bw[mindex2(1,z,4)] = ( 3*pow(w,3) - 6*pow(w,2) + 4)/6;
    Bw[mindex2(2,z,4)] = (-3*pow(w,3) + 3*pow(w,2) + 3*w + 1)/6;
    Bw[mindex2(3,z,4)] = pow(w,3)/6;
  }
  
  ZOsize_d[0]=ZOsizex;  ZOsize_d[1]=ZOsizey; ZOsize_d[2]=ZOsizez;
  Osize_d[0]=Osizex;  Osize_d[1]=Osizey; Osize_d[2]=Osizez;
  
  nlhs_d[0]=(double)nlhs;
  
  WaitForThread1 = 1;
  ThreadArgs1[0]=Bu;
  ThreadArgs1[1]=Bv;
  ThreadArgs1[2]=Bw;
  ThreadArgs1[3]=ZOsize_d;
  ThreadArgs1[4]=Osize_d;
  ThreadArgs1[5]=ZI;
  ThreadArgs1[6]=Tx;
  ThreadArgs1[7]=Ty;
  ThreadArgs1[8]=Tz;
  ThreadArgs1[9]=dxa;
  ThreadArgs1[10]=dya;
  ThreadArgs1[11]=dza;
  ThreadArgs1[12]=ThreadID1;
  ThreadArgs1[13]=Ox;
  ThreadArgs1[14]=Oy;
  ThreadArgs1[15]=Oz;
  ThreadArgs1[16]=ZO;
  ThreadArgs1[17]=nlhs_d;
  ThreadList[0] = (HANDLE)_beginthreadex( NULL, 0, &transformvolume, ThreadArgs1 , 0, NULL );
    
  WaitForThread2 = 1;
  ThreadArgs2[0]=Bu;
  ThreadArgs2[1]=Bv;
  ThreadArgs2[2]=Bw;
  ThreadArgs2[3]=ZOsize_d;
  ThreadArgs2[4]=Osize_d;
  ThreadArgs2[5]=ZI;
  ThreadArgs2[6]=Tx;
  ThreadArgs2[7]=Ty;
  ThreadArgs2[8]=Tz;
  ThreadArgs2[9]=dxa;
  ThreadArgs2[10]=dya;
  ThreadArgs2[11]=dza;
  ThreadArgs2[12]=ThreadID2;
  ThreadArgs2[13]=Ox;
  ThreadArgs2[14]=Oy;
  ThreadArgs2[15]=Oz;
  ThreadArgs2[16]=ZO;
  ThreadArgs2[17]=nlhs_d;
  ThreadList[1] = (HANDLE)_beginthreadex( NULL, 0, &transformvolume, ThreadArgs2 , 0, NULL );

  
  WaitForSingleObject(ThreadList[0], INFINITE);
  WaitForSingleObject(ThreadList[1], INFINITE);
  
  CloseHandle( ThreadList[0] );
  CloseHandle( ThreadList[1] ); 
}
        

