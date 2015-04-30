/*
* Yongxiang Huang, last modification: 28/03/2010
* yongxianghuang@gmail.com
*
*/

/* This function is to estimate the structure function scaling by considering a sign power
 * The  structure function is further termed into four part
 * S_1(\tau) (sfP) is the structure function with positive contribution
 * S_2(\tau) (sfN) is the structure function with negative contribution
 * S_3(\tau) (sfO) is the structure function with all contribution (|P|+|N|)
 * S_4(\tau) (sfM) is the structure function with all contribution (|P|-|N|)
 */

#include <stdlib.h>
#include <stdio.h>
#include "mex.h"
#include <math.h>






/************************************************************************/
/*                                                                      */
/* MAIN FUNCTION                                                        */
/*                                                                      */
/************************************************************************/

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
  
    /* declarations */
    int i,j,m,n,k,k1,k2,k3,k4,k5,tn,Nx,Ntau,Nq;
    
    double *dq,*sfN,*sfP,*sfO,*sfM,*tau,*tmp,*tmpT,temp,*Nk,*q,*Xsign,*x; 
    /*instantaneous frequency, instantaneous amplitude, resolution of amplitude resolution frequency*/
      
    
/*     check input*/
    if (nrhs!=4)     mexErrMsgTxt("You have to input two parameters!");
    if (mxIsEmpty(prhs[0]))mexErrMsgTxt("Time series is empty!");
    if (mxIsEmpty(prhs[1]))mexErrMsgTxt("The time delay is empty!");
    if (mxIsEmpty(prhs[2]))mexErrMsgTxt("q is empty!");
    if (mxIsEmpty(prhs[3]))mexErrMsgTxt("dq is empty!");
    /* get input data */
    x=mxGetPr(prhs[0]);
    tau=mxGetPr(prhs[1]);
    q=mxGetPr(prhs[2]);
    dq=mxGetPr(prhs[3]);
            
    Nx=mxGetN(prhs[0]);
    tn=mxGetM(prhs[0]);
    if (tn>Nx) Nx=tn;
    
    Ntau=mxGetN(prhs[1]);/*length of time delay*/
    tn=mxGetM(prhs[1]);
    if (tn>Ntau) Ntau=tn;
    Nq=q[0]/dq[0];
    Xsign=(double *)malloc(Nx*sizeof(double)); /*specify sign*/
    tmp=(double *)malloc(Nx*sizeof(double)); /*specify tmp for moments*/
    tmpT=(double *)malloc(Nx*sizeof(double)); /*specify tmp for moments*/
    
    plhs[0]=mxCreateDoubleMatrix(Nq,Ntau,mxREAL);
    sfP=mxGetPr(plhs[0]);
    plhs[1]=mxCreateDoubleMatrix(Nq,Ntau,mxREAL);
    sfN=mxGetPr(plhs[1]);
    plhs[2]=mxCreateDoubleMatrix(Nq, Ntau, mxREAL);
    sfO=mxGetPr(plhs[2]);
    plhs[3]=mxCreateDoubleMatrix(Nq,Ntau,mxREAL);
    sfM=mxGetPr(plhs[3]);
    
    plhs[4]=mxCreateDoubleMatrix(2,Ntau,mxREAL);
    Nk=mxGetPr(plhs[4]);
     
  
  
    for(i=0;i<Ntau;i++)/*Bigest loop for for different time delay*/
    {
        k2=tau[i];  /*separation tau*/
        k1=Nx-k2; /*length of velocity average*/
      
         
        for(j=0;j<k1;j++) /*velocity average*/
        {
            tmp[j]=(x[j+k2]+x[j])/2;
           /* if(tmp[j]>0){ Xsign[j]=1.0;}
            if(tmp[j]<0){ Xsign[j]=-1.0;tmp[j]=-tmp[j];}
            if(tmp[j]==0){ Xsign[j]=0.0;}
            tmp[j]=pow(tmp[j], dq[0]);
            tmpT[j]=tmp[j];*/
        }
        
     /*   temp=0.0;/*remove mean value*/
     /*   for(j=0;j<k1;j++)temp=temp+tmp[j] ;
        temp=temp/k1;
        for(j=0;j<k1;j++)tmp[j]=tmp[j]-temp;*/
        
        for(j=0;j<k1;j++) /*velocity sign and sign power*/
        {
            if(tmp[j]>0){ Xsign[j]=1.0;}
            if(tmp[j]<0){ Xsign[j]=-1.0;tmp[j]=-tmp[j];}
            if(tmp[j]==0){ Xsign[j]=0.0;}
            tmp[j]=pow(tmp[j], dq[0]);
            tmpT[j]=tmp[j];
        }
        
        k4=0;
        k5=0;
       
        for(k=0;k<Nq;k++) /*loop for different q*/
        {
            for(j=0;j<k1;j++)/*loop for autocorrelation function*/
            {
                
                if(Xsign[j]>0)
                {
                    sfP[i*Nq+k]=sfP[i*Nq+k]+tmpT[j];
                    k4=k4+1;
                }
                if(Xsign[j]<0)
                {
                    sfN[i*Nq+k]=sfN[i*Nq+k]+tmpT[j];
                    k5=k5+1;
                }
                sfO[i*Nq+k]=sfO[i*Nq+k]+tmpT[j];
            }/*loop for autocorrelation function*/
           for(j=0;j<k1;j++)tmpT[j]=tmpT[j]*tmp[j];
           
        }/*loop for different q*/
        Nk[i*2]=k4/Nq;
        Nk[i*2+1]=k5/Nq;
       /*    mexPrintf("OK here tau[%i] = %i \n", i+1,k2);
       if(i==2)return;*/
        for(k=0;k<Nq;k++)
        {
            sfM[i*Nq+k]=(sfP[i*Nq+k]-sfN[i*Nq+k])/k1;
            if(Nk[i*2]>0)sfP[i*Nq+k]=sfP[i*Nq+k]/Nk[i*2];
            if(Nk[i*2+1]>0)sfN[i*Nq+k]=sfN[i*Nq+k]/Nk[i*2+1];
            sfO[i*Nq+k]=sfO[i*Nq+k]/k1;
            
        }
       
      
  }
    mxFree(Xsign);
    mxFree(tmp);
    mxFree(tmpT);
}