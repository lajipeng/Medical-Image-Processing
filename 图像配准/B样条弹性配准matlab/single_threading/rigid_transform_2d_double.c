#include "mex.h"
#include "math.h"

/* Imrigide transform Iout=rigide_transform_2d_double(Iin,Minv)
  * This function transforms a image with a 3x3 transformation matrix 
  * 
  * Function is written by D.Kroon University of Twente (September 2008)
  */

// Convert 2D/3D matrix index to 1D index
int mindex2(int x, int y, int sizx) { return y*sizx+x; }

// The matlab mex function
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    // Ox and Oy are the grid points
    // Zo is the input image
    // Zi is the transformed image
    // nx and ny are the number of grid points (inside the image)
    double *Iin, *Iout, *M;
    
    
    // Transformation matrix
    double A00, A10, A20;
    double A01, A11, A21;
    double A02, A12, A22;

    // Size of input image
    mwSize  Isizex, Isizey;
    const mwSize *dims;
    
    // Variables to store 1D index
    int indexI;
    
    // Location of pixel which will be come the current pixel
    double Tlocalx;
    double Tlocaly;
    
    // Linear interpolation variables
    int xBas[4], yBas[4];
    double perc[4];
    double xCom, yCom;
    double color[4]={0,0,0,0};
    
    // X,Y coordinates of current pixel
    int x,y;
    double xmean, ymean;
    double xd,yd;
    
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

  // Create output array
  plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);

  /* Assign pointers to each input. */
  Iin=mxGetPr(prhs[0]);
  M=mxGetPr(prhs[1]);
 
  A00 = M[mindex2(0,0,3)]; A10 = M[mindex2(0,1,3)]; A20 = M[mindex2(0,2,3)]; 
  A01 = M[mindex2(1,0,3)]; A11 = M[mindex2(1,1,3)]; A21 = M[mindex2(1,2,3)]; 
  A02 = M[mindex2(2,0,3)]; A12 = M[mindex2(2,1,3)]; A22 = M[mindex2(2,2,3)]; 
  
  /* Assign pointer to output. */
  Iout = mxGetPr(plhs[0]);
    
  xmean=((double)Isizex)/2;
  ymean=((double)Isizey)/2;
  
  // Loop through all image pixel coordinates
  for (y=0; y<Isizey; y++)
  {
    for (x=0; x<Isizex; x++)
    {
        xd=(double)x-xmean; yd=(double)y-ymean;
                 
        Tlocalx = xmean + A00 * xd + A10 *yd + A20 * 1;
        Tlocaly = ymean + A01 * xd + A11 *yd + A21 * 1;
                                
        // Determine the coordinates of the pixel(s) which will be come the current pixel
        // (using linear interpolation)  
        xBas[0]=(int) floor(Tlocalx); yBas[0]=(int) floor(Tlocaly);
        xBas[1]=xBas[0]+0;      yBas[1]=yBas[0]+1;
        xBas[2]=xBas[0]+1;      yBas[2]=yBas[0]+0;
        xBas[3]=xBas[0]+1;      yBas[3]=yBas[0]+1;

        if(xBas[0]>=0&&xBas[0]<Isizex&&yBas[0]>=0&&yBas[0]<Isizey) { indexI=mindex2(xBas[0],yBas[0],Isizex); color[0]=Iin[indexI]; } else { color[0]=0;}
        if(xBas[1]>=0&&xBas[1]<Isizex&&yBas[1]>=0&&yBas[1]<Isizey) { indexI=mindex2(xBas[1],yBas[1],Isizex); color[1]=Iin[indexI]; } else { color[1]=0;}
        if(xBas[2]>=0&&xBas[2]<Isizex&&yBas[2]>=0&&yBas[2]<Isizey) { indexI=mindex2(xBas[2],yBas[2],Isizex); color[2]=Iin[indexI]; } else { color[2]=0;}
        if(xBas[3]>=0&&xBas[3]<Isizex&&yBas[3]>=0&&yBas[3]<Isizey) { indexI=mindex2(xBas[3],yBas[3],Isizex); color[3]=Iin[indexI]; } else { color[3]=0;}

        // Linear interpolation constants (percentages)
        xCom=Tlocalx-floor(Tlocalx); yCom=Tlocaly-floor(Tlocaly);
        perc[0]=(1-xCom) * (1-yCom);
        perc[1]=(1-xCom) * yCom;
        perc[2]=xCom * (1-yCom);
        perc[3]=xCom * yCom;

        // Set the current pixel value
        indexI=mindex2(x,y,Isizex);
        Iout[indexI]=color[0]*perc[0]+color[1]*perc[1]+color[2]*perc[2]+color[3]*perc[3];
    }
  }
}
        

