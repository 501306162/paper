#include "mex.h"
#include "math.h"

/* This function makes a 2D joint histogram of 1D,2D...ND images
 * and also calculates the seperate histograms of both images.
 *
 * [hist12, hist1, hist2]=mutual_histogram_double(I1,I2,Imin,Imax,nbins);
 *
 * Function is written by D.Kroon University of Twente (July 2008)
 */

int mindex2(int x, int y, int sizx) { return y*sizx+x; }

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    /* I1 and I2 are the input images */
    /* hist12 joint histogram */
    /* hist1 and histogram of I1 and hist2 of I2 */
    double *I1, *I2, *Imin, *Imax, *nbins, *hist12, *hist1, *hist2;
    
    /* Size of input image */
    const mwSize *idims; 
    
    /* Size of output image */
    int odims2[2]={0,0};
    int odims1[1]={0};    
    /* Dimensions */
    int nsubs;
    int npixels=1;
    
    /* index var */
    int index;
    
    /* intensity location */
    int sizex;
    double xd, xm, xp, xmd, xpd;
    double yd, ym, yp, ymd, ypd;
        
    
    /* loop vars */
    int i;
    
    /* vars */
    double minv;
    double scav;
    
    /* Check for proper number of arguments. */
    if(nrhs!=5) {
       mexErrMsgTxt("five inputs are required.");
    } else if(nlhs!=3) {
       mexErrMsgTxt("Three outputs are required");
    }
  
    /* Get the number of dimensions */
    nsubs = mxGetNumberOfDimensions(prhs[0]);
    /* Get the sizes of the grid */
    idims = mxGetDimensions(prhs[0]);   
    for (i=0; i<nsubs; i++) { npixels=npixels*idims[i]; }

    /* Assign pointers to each input. */
    I1=(double *)mxGetData(prhs[0]);
    I2=(double *)mxGetData(prhs[1]);
    Imin=(double *)mxGetData(prhs[2]);
    Imax=(double *)mxGetData(prhs[3]);
    nbins=(double *)mxGetData(prhs[4]);
    
    /* Create image matrix for the return arguments */
    odims2[0]=(int) nbins[0]; odims2[1]=(int)nbins[0];  
    plhs[0] = mxCreateNumericArray(2, odims2, mxDOUBLE_CLASS, mxREAL);
    odims1[0]=(int) nbins[0]; 
    plhs[1] = mxCreateNumericArray(1, odims1, mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(1, odims1, mxDOUBLE_CLASS, mxREAL);

    /* Assign pointers to each output. */
    hist12=(double *)mxGetData(plhs[0]);
    hist1=(double *)mxGetData(plhs[1]);
    hist2=(double *)mxGetData(plhs[2]);

    /* min value */
    minv=Imin[0];
    /* scale value */
    scav=nbins[0]/(Imax[0]-Imin[0]);
    sizex=(int) nbins[0];
    for (i=0; i<npixels; i++)
    {
        xd=(double)scav*(I1[i]-minv);
        xm=floor(xd); xp=xm+1;
        xmd=xp-xd; xpd=xd-xm;
                
        yd=(double)scav*(I2[i]-minv);
        ym=floor(yd); yp=ym+1;
        ymd=yp-yd; ypd=yd-ym;
        
        if(xm<0){ xm=0; } else if(xm>(double)(sizex-1)) { xm=(double)(sizex-1); }
        if(xp<0){ xp=0; } else if(xp>(double)(sizex-1)) { xp=(double)(sizex-1); }
        if(ym<0){ ym=0; } else if(ym>(double)(sizex-1)) { ym=(double)(sizex-1); }
        if(yp<0){ yp=0; } else if(yp>(double)(sizex-1)) { yp=(double)(sizex-1); }
        
        index=mindex2((int)xm,(int)ym,sizex);
        hist12[index]=hist12[index]+xmd*ymd;
        index=mindex2((int)xp,(int)ym,sizex);
        hist12[index]=hist12[index]+xpd*ymd;
        index=mindex2((int)xm,(int)yp,sizex);
        hist12[index]=hist12[index]+xmd*ypd;
        index=mindex2((int)xp,(int)yp,sizex);
        hist12[index]=hist12[index]+xpd*ypd;

        hist1[(int)xm]=hist1[(int)xm]+xmd; hist1[(int)xp]=hist1[(int)xp]+xpd;
        hist2[(int)ym]=hist2[(int)ym]+ymd; hist2[(int)yp]=hist2[(int)yp]+ypd;
    }
}
        

