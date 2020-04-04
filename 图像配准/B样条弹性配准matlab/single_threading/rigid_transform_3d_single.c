#include "mex.h"
#include "math.h"

/* Imrigide transform Iout=rigide_transform_3d_single(Iin,Minv)
  * This function transforms a volume with a 4x4 transformation matrix 
  * 
  * Function is written by D.Kroon University of Twente (September 2008)
  */

// Convert 2D/3D matrix index to 1D index
int mindex2(int x, int y, int sizx) { return y*sizx+x; }
int mindex3(int x, int y, int z, int sizx, int sizy) { return z*sizx*sizy+y*sizx+x;}

// The matlab mex function
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    // Ox and Oy are the grid points
    // Zo is the input image
    // Zi is the transformed image
    // nx and ny are the number of grid points (inside the image)
    float *Iin, *Iout, *M;
    
    
    // Transformation matrix
    float A00, A10, A20, A30;
    float A01, A11, A21, A31;
    float A02, A12, A22, A32;
    float A03, A13, A23, A33;

    // Size of input image
    mwSize  Isizex, Isizey, Isizez;
    const mwSize *dims;
    
    // Variables to store 1D index
    int indexI;
    
    // loop variables
    int p;
    
    // Location of pixel which will be come the current pixel
    float Tlocalx;
    float Tlocaly;
    float Tlocalz;
    
    // Linear interpolation variables
    int xBas[8], yBas[8], zBas[8];
    float perc[8];
    float xCom, yCom, zCom;
    float color[8]={0,0,0,0,0,0,0,0};
    
    // X,Y,Z coordinates of current pixel
    int x,y,z;
    float xmean, ymean, zmean;
    float xd,yd,zd;
    
  /* Check for proper number of arguments. */
  if(nrhs!=2) {
    mexErrMsgTxt("Two inputs are required.");
  } else if(nlhs!=1) {
    mexErrMsgTxt("One output required");
  }
    
  // Get the sizes of the image
  dims = mxGetDimensions(prhs[0]);   
  Isizex = dims[0]; 
  Isizey = dims[1];
  Isizez = dims[2];
  
  // Create output array
  plhs[0] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);

  /* Assign pointers to each input. */
  Iin=(float *)mxGetData(prhs[0]);
  M=(float *)mxGetData(prhs[1]);
 
  A00 = M[mindex2(0,0,4)]; A10 = M[mindex2(0,1,4)]; A20 = M[mindex2(0,2,4)]; A30 = M[mindex2(0,3,4)];
  A01 = M[mindex2(1,0,4)]; A11 = M[mindex2(1,1,4)]; A21 = M[mindex2(1,2,4)]; A31 = M[mindex2(1,3,4)];
  A02 = M[mindex2(2,0,4)]; A12 = M[mindex2(2,1,4)]; A22 = M[mindex2(2,2,4)]; A32 = M[mindex2(2,3,4)];
  A03 = M[mindex2(3,0,4)]; A13 = M[mindex2(3,1,4)]; A23 = M[mindex2(3,2,4)]; A33 = M[mindex2(3,3,4)];
  
  /* Assign pointer to output. */
  Iout = (float *)mxGetData(plhs[0]);
    
  xmean=((float)Isizex)/2;
  ymean=((float)Isizey)/2;
  zmean=((float)Isizez)/2;
  
  // Loop through all image pixel coordinates
      for (z=0; z<Isizez; z++)
      {
        for (y=0; y<Isizey; y++)
        {
            for (x=0; x<Isizex; x++)
            {
                xd=(float)x-xmean; yd=(float)y-ymean; zd=(float)z-zmean;
                
                Tlocalx = xmean + A00 * xd + A10 *yd + A20 *zd + A30 * 1;
                Tlocaly = ymean + A01 * xd + A11 *yd + A21 *zd + A31 * 1;
                Tlocalz = zmean + A02 * xd + A12 *yd + A22 *zd + A32 * 1;
                                
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
                    if(xBas[p]>=0&&xBas[p]<Isizex&&yBas[p]>=0&&yBas[p]<Isizey&&zBas[p]>=0&&zBas[p]<Isizez) 
                    { 
                        indexI=mindex3(xBas[p],yBas[p],zBas[p],Isizex,Isizey); 
                        color[p]=Iin[indexI]; 
                    } 
                    else { color[p]=0; }
                }

                // Linear interpolation constants (percentages)
                xCom=Tlocalx-(float)floor((double)Tlocalx); yCom=Tlocaly-(float)floor((double)Tlocaly);  zCom=Tlocalz-(float)floor((double)Tlocalz);
                perc[0]=(1-xCom) * (1-yCom) * (1-zCom);
                perc[1]=(1-xCom) * (1-yCom) * zCom;
                perc[2]=(1-xCom) * yCom * (1-zCom);
                perc[3]=(1-xCom) * yCom * zCom;
                perc[4]=xCom * (1-yCom) * (1-zCom);
                perc[5]=xCom * (1-yCom) * zCom;
                perc[6]=xCom * yCom * (1-zCom);
                perc[7]=xCom * yCom * zCom;

                // Set the current pixel value
                indexI=mindex3(x,y,z,Isizex,Isizey);
                Iout[indexI]=0;
                for (p=0; p<8; p++)
                {
                    Iout[indexI]+=color[p]*perc[p];
                }
            }
        }
    }
}
        

