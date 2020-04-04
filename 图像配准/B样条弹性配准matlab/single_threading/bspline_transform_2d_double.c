#include "mex.h"
#include "math.h"

/* Bspline transformation grid function
 * function [Iout,Tx,Ty]=bspline_transform_2d_double(Ox,Oy,Iin,dx,dy)
 * 
 * Ox, Oy are the grid points coordinates
 * Iin is input image, Iout the transformed output image
 * dx and dy are the spacing of the b-spline knots
 *
 * Iout: The transformed image
 * Tx: The transformation field in x direction
 * Ty: The transformation field in y direction
 *
 * This function is an implementation of the b-spline registration
 * algorithm in "D. Rueckert et al. : Nonrigid Registration Using Free-Form 
 * Deformations: Application to Breast MR Images".
 * 
 * We used "Fumihiko Ino et al. : a data distrubted parallel algortihm for 
 * nonrigid image registration" for the correct formula's, because 
 * (most) other papers contain errors. 
 *
 *   Function is written by D.Kroon University of Twente (July 2008)
 */

// Convert 2D matrix index to 1D index
int mindex2(int x, int y, int sizx, int sizy) 
{
    if(x<0) { x=0; }
    if(x>(sizx-1)) { x=sizx-1; }
    if(y<0) { y=0; }
    if(y>(sizy-1)) { y=sizy-1; }
    return y*sizx+x;
}

// The matlab mex function
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    // Ox and Oy are the grid points
    // Zo is the input image
    // Zi is the transformed image
    // dx and dy are the spacing of the b-spline knots
    double *Ox,*Oy, *ZO, *dxa, *dya, *ZI,*Tx,*Ty;
    
    // Size of input image
    mwSize  ZOsizex, ZOsizey;
    // Size of grid
    mwSize  Osizex, Osizey;
    
    // B-spline variables
    double u,v;
    int u_index=0; 
    int v_index=0;
    int i, j;
    double *Bu, *Bv;
    
    // B-Spline loop variabels
    int l,m;
   
    // Grid distance;
    int dx,dy; 
    
    // Location of pixel which will be come the current pixel
    double Tlocalx;
    double Tlocaly;
    
    // Variables to store 1D index
    int indexO;
    int indexZI;
    int indexZO;
     
    // X,Y coordinates of current pixel
    int x,y;
    
    // Linear interpolation variables
    int xBas[4], yBas[4];
    double perc[4];
    double xCom, yCom;
    double color[4]={0,0,0,0};
             
  
  /* Check for proper number of arguments. */
  if(nrhs!=5) {
    mexErrMsgTxt("Five inputs are required.");
  }
   
  // Get the sizes of the grid
  Osizex = (mwSize)mxGetM(prhs[0]);  
  Osizey = (mwSize)mxGetN(prhs[0]);
  
  // Create image matrix for the return arguments with the size of input image   
  ZOsizex = (mwSize) mxGetM(prhs[2]);  ZOsizey = (mwSize) mxGetN(prhs[2]);
  plhs[0] = mxCreateDoubleMatrix(ZOsizex,ZOsizey, mxREAL); 
  if(nlhs>1) { plhs[1] = mxCreateDoubleMatrix(ZOsizex,ZOsizey, mxREAL); }
  if(nlhs>2) { plhs[2] = mxCreateDoubleMatrix(ZOsizex,ZOsizey, mxREAL); }
  
  
  /* Assign pointers to each input. */
  Ox=mxGetPr(prhs[0]);
  Oy=mxGetPr(prhs[1]);
  ZO=mxGetPr(prhs[2]);
  dxa=mxGetPr(prhs[3]);
  dya=mxGetPr(prhs[4]);
  
  /* Get the spacing of the uniform b-spline grid */
  dx=(int)dxa[0]; dy=(int)dya[0];
  
  
  /* Assign pointer to output. */
  ZI = mxGetPr(plhs[0]);
  if(nlhs>1) { Tx = mxGetPr(plhs[1]); }
  if(nlhs>2) { Ty = mxGetPr(plhs[2]); }
  
  // Make polynomial look up tables 
  Bu=malloc(dx*4*sizeof(double));
  Bv=malloc(dy*4*sizeof(double));
  for (x=0; x<dx; x++)
  {
    u=(x/(double)dx)-floor(x/(double)dx);
    Bu[mindex2(0,x,4,dx)] = pow((1-u),3)/6;
    Bu[mindex2(1,x,4,dx)] = ( 3*pow(u,3) - 6*pow(u,2) + 4)/6;
    Bu[mindex2(2,x,4,dx)] = (-3*pow(u,3) + 3*pow(u,2) + 3*u + 1)/6;
    Bu[mindex2(3,x,4,dx)] = pow(u,3)/6;
  }
  
  for (y=0; y<dy; y++)
  {
    v=(y/(double)dy)-floor(y/(double)dy);
    Bv[mindex2(0,y,4,dy)] = pow((1-v),3)/6;
    Bv[mindex2(1,y,4,dy)] = ( 3*pow(v,3) - 6*pow(v,2) + 4)/6;
    Bv[mindex2(2,y,4,dy)] = (-3*pow(v,3) + 3*pow(v,2) + 3*v + 1)/6;
    Bv[mindex2(3,y,4,dy)] = pow(v,3)/6;
  }
  
  
    
  // Loop through all image pixel coordinates
    for (y=0; y<ZOsizey; y++)
    {
        for (x=0; x<ZOsizex; x++)
        {
            // Calculate the indexes need to loop up the B-spline values.
            u_index=x%dx; 
            v_index=y%dy;
            
            i=(int)floor(x/(double)dx); // (first row outside image against boundary artefacts)
            j=(int)floor(y/(double)dy);
        
            // This part calculates the coordinates of the pixel
            // which will be transformed to the current x,y pixel.
            Tlocalx=0; Tlocaly=0;
            for(l=0; l<4; l++)
            {
                for(m=0; m<4; m++)
                {    
                     if(((i+l)>=0)&&((i+l)<Osizex)&&((j+m)>=0)&&((j+m)<Osizey))
                     {
                          indexO=mindex2(i+l,j+m,Osizex,Osizey);
                          Tlocalx=Tlocalx+Bu[mindex2(l,u_index,4,dx)]*Bv[mindex2(m,v_index,4,dy)]*Ox[indexO];
                          Tlocaly=Tlocaly+Bu[mindex2(l,u_index,4,dx)]*Bv[mindex2(m,v_index,4,dy)]*Oy[indexO];
                     }
                }
            }            

            // Determine the coordinates of the pixel(s) which will be come the current pixel
            // (using linear interpolation)  
            xBas[0]=(int) floor(Tlocalx); yBas[0]=(int) floor(Tlocaly);
            xBas[1]=xBas[0]+0;      yBas[1]=yBas[0]+1;
            xBas[2]=xBas[0]+1;      yBas[2]=yBas[0]+0;
            xBas[3]=xBas[0]+1;      yBas[3]=yBas[0]+1;

            if(xBas[0]>=0&&xBas[0]<ZOsizex&&yBas[0]>=0&&yBas[0]<ZOsizey) { indexZI=mindex2(xBas[0],yBas[0],ZOsizex,ZOsizey); color[0]=ZO[indexZI]; } else { color[0]=0;}
            if(xBas[1]>=0&&xBas[1]<ZOsizex&&yBas[1]>=0&&yBas[1]<ZOsizey) { indexZI=mindex2(xBas[1],yBas[1],ZOsizex,ZOsizey); color[1]=ZO[indexZI]; } else { color[1]=0;}
            if(xBas[2]>=0&&xBas[2]<ZOsizex&&yBas[2]>=0&&yBas[2]<ZOsizey) { indexZI=mindex2(xBas[2],yBas[2],ZOsizex,ZOsizey); color[2]=ZO[indexZI]; } else { color[2]=0;}
            if(xBas[3]>=0&&xBas[3]<ZOsizex&&yBas[3]>=0&&yBas[3]<ZOsizey) { indexZI=mindex2(xBas[3],yBas[3],ZOsizex,ZOsizey); color[3]=ZO[indexZI]; } else { color[3]=0;}
  
            // Linear interpolation constants (percentages)
            xCom=Tlocalx-floor(Tlocalx); yCom=Tlocaly-floor(Tlocaly);
            perc[0]=(1-xCom) * (1-yCom);
            perc[1]=(1-xCom) * yCom;
            perc[2]=xCom * (1-yCom);
            perc[3]=xCom * yCom;

            // Set the current pixel value
            indexZO=mindex2(x,y,ZOsizex,ZOsizey);
            ZI[indexZO]=color[0]*perc[0]+color[1]*perc[1]+color[2]*perc[2]+color[3]*perc[3];
            
            // Store transformation field
            if(nlhs>1) { Tx[indexZO]=Tlocalx-(double)x; }
            if(nlhs>2) { Ty[indexZO]=Tlocaly-(double)y; }
        }
    }
}
        

