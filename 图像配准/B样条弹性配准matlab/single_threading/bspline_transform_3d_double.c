#include "mex.h"
#include "math.h"

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
    double *Ox,*Oy,*Oz,*ZO, *dxa, *dya,*dza,*ZI, *Tx,*Ty,*Tz;
    
    // Size of input image
    mwSize  ZOsizex, ZOsizey, ZOsizez;
    const mwSize *dims;
    // Size of grid
    mwSize  Osizex, Osizey, Osizez;
    
    // B-spline variables
    double u,v,w;
    int u_index=0; 
    int v_index=0;
    int w_index=0;
    
    int i,j,k;
    double *Bu, *Bv, *Bw;
    
    // B-Spline loop variabels
    int l,m,n;
   
    // loop variable
    int p;
    
    // Grid distance;
    int dx,dy,dz; 
    
    // Location of pixel which will be come the current pixel
    double Tlocalx;
    double Tlocaly;
    double Tlocalz;
        
    // Variables to store 1D index
    int indexO;
    int indexZI;
    int indexZO;
     
    // X,Y,Z coordinates of current pixel
    int x,y,z;
    
    // Linear interpolation variables
    int xBas[8], yBas[8], zBas[8];
    double perc[8];
    double xCom, yCom, zCom;
    double color[8]={0,0,0,0,0,0,0,0};
             
  
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
  
  // Loop through all image pixel coordinates
     for (z=0; z<ZOsizez; z++)
     {
        for (y=0; y<ZOsizey; y++)
        {
            for (x=0; x<ZOsizex; x++)
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
                            if(((i+l)>=0)&&((i+l)<Osizex)&&((j+m)>=0)&&((j+m)<Osizey)&&((k+n)>=0)&&((k+n)<Osizez))
                            {
                                 indexO=mindex3(i+l,j+m,k+n,Osizex,Osizey);
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
                    if(xBas[p]>=0&&xBas[p]<ZOsizex&&yBas[p]>=0&&yBas[p]<ZOsizey&&zBas[p]>=0&&zBas[p]<ZOsizez) 
                    { 
                        indexZI=mindex3(xBas[p],yBas[p],zBas[p],ZOsizex,ZOsizey); 
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
                indexZO=mindex3(x,y,z,ZOsizex,ZOsizey);
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
}
        

