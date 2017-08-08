#include <stdio.h>

__device__ void alpha_func(double *p, double dt_chem, double *alpha, int nspec, int idx, int numcell)               
{                                                          
    double rr;                                                 
    for(int i=0; i < nspec; i++) {                                                          
        rr = 1.0/(p[idx+numcell*i]*dt_chem+1.0e-50);        
        alpha[idx+numcell*i] = (180.0*rr*rr*rr+60.0*rr*rr+11.0*rr+1.0)/(360.0*rr*rr*rr+60.0*rr*rr+12.0*rr+1.0);     
    }                                                          
}                                                          

__device__ void ydot(double *Y, double *MM, double *kf, double *kb, double *p, double *q, int idx, int numcell)       
{            
// <insert-source-terms-here>
}                                                          

__global__ void new_qss(double *kf, double *kb, double *MM, double *Y, double *Yp, double *Yc,
			double *h, double *alpha, double *q, double *p, double *qp, double *pp,
			double *q_bar, double *p_bar, double *alpha_bar, double *Ypo, const int nspec,
			const  int numreac, const  int numcell, double dt_global, double *debugging)                   
{                                                          
double dt_flow = dt_global;                                     
double t = 0.0;                                                  
int flag;                                                                   
double e = 0.001;                                                           
double epsilon = e/2.0;                                                      
int step_fail;                                                             
int idx = blockIdx.x;                                                        
double dt_chem;                                                              
dt_chem = h[idx];                                                                     
int max_iterations = 10;                                                     
int count = 0;                                                               
double x0 = 0;                                                            
double cond = 0.0;                                                         
double sigma = 0.0;                                                        
double dt_chem_temp = 0;
while (t < dt_flow) {
    count = 0;                                                             
    flag = 0;                                                                       
    step_fail = 1;                                                       
    ydot(Y, MM, kf, kb, p, q, idx, numcell);                                        
    alpha_func(p, dt_chem, alpha, nspec, idx, numcell);             
    for (int i = 0; i < nspec; i++) {                                                                        
        Yp[idx+numcell*i] = Y[idx+numcell*i] +
	  (dt_chem*(q[idx+numcell*i]-p[idx+numcell*i]*Y[idx+numcell*i]))/(1.0+alpha[idx+numcell*i]*dt_chem*p[idx+numcell*i]);
        Ypo[idx+numcell*i] = Yp[idx+numcell*i];      
    }                                                                        
    while (step_fail == 1) {                                                                          
        step_fail = 0;  
        sigma = 1.0e-50;   
        cond = 0.0; 
        ydot(Yp, MM, kf, kb, pp, qp, idx, numcell);      
        for (int j = 0; j < nspec; j++) {            
            p_bar[idx+numcell*j] = 0.5*(pp[idx+numcell*j]+p[idx+numcell*j]); 
        }  
        alpha_func(p_bar, dt_chem, alpha_bar, nspec, idx, numcell); 
        for (int r=0; r < nspec; r++) {    
            q_bar[idx+numcell*r] = alpha_bar[idx+numcell*r]*qp[idx+numcell*r]+(1.0-alpha_bar[idx+numcell*r])*q[idx+numcell*r]; 
        }     
        for (int k = 0; k < nspec; k++) { 
            Yc[idx+numcell*k] = Y[idx+numcell*k] +
	      (dt_chem*(q_bar[idx+numcell*k]-p_bar[idx+numcell*k]*Y[idx+numcell*k]))/
	      (1.0+alpha_bar[idx+numcell*k]*dt_chem*p_bar[idx+numcell*k]);
        }    
        for (int l = 0; l < nspec; l++) { 
            cond = (fabs(Yc[idx+numcell*l]-Ypo[idx+numcell*l])); 
            if (cond >= (e*(Yc[idx+numcell*l]+1.0e-10))) {        
                step_fail = 1;      
            }       
        }     
        count = count + 1;    
        if (count > max_iterations) {   
            step_fail = 0; 
            flag  = 1; 
        } 
        if (step_fail == 1) {         
            for (int u = 0; u < nspec; u++) {
	        Yp[idx+numcell*u] = Yc[idx+numcell*u];       
            }      
        }   
    }
    if (flag == 0) {   
        for (int m = 0; m < nspec; m++) {
            Y[idx+numcell*m] = Yc[idx+numcell*m];
        }  
        t = t + dt_chem; 
        h[idx] = dt_chem_temp; 
    } 	
    for (int n=0; n < nspec; n++) { 
        cond = (fabs(Yc[idx+numcell*n]-Ypo[idx+numcell*n])/(epsilon*Yc[idx+numcell*n])); 
        if (cond > sigma) {    
            sigma = cond;
        } 
    }
    if (sigma <= 0.0) {
        dt_chem_temp = dt_chem;
    }
    else {
        x0 = sigma;
        x0 = x0 - 0.5*(x0*x0 - sigma)/x0;
        x0 = x0 - 0.5*(x0*x0 - sigma)/x0;
        x0 = x0 - 0.5*(x0*x0 - sigma)/x0;
        dt_chem_temp = dt_chem * ((1.0/x0) + 0.005);
    }
    if (flag == 0) {   
        if(dt_chem_temp/dt_chem > 1.15) {     
            dt_chem = 1.15*dt_chem;   
        }    
        else {     
            dt_chem = dt_chem_temp;        
        }       
    }       
    if (flag == 1) {    
        if(dt_chem_temp/dt_chem > (1.0/3.0)) {        
            dt_chem = (1.0/3.0)*dt_chem;    
        }     
        else if(dt_chem_temp/dt_chem < 0.01) {    
            dt_chem = 0.01*dt_chem;
        }  
        else {
            dt_chem = dt_chem_temp;
        }
    }
    dt_chem_temp = dt_chem;  
    dt_chem = fmin(dt_flow - t, dt_chem); 
    }    
} 


extern "C" {
       void kernel_launcher(double *kf, double *kb, double *MM, double *Y, double *Yp, double *Yc, 
       	    			   double *h, double *alpha, double *q, double *p, double *qp, double *pp, 
				   double *q_bar, double *p_bar, double *alpha_bar, double *Ypo, const int nspec, 
				   const  int numreac, const  int numcell, double dt_global, double *debugging, 
				   size_t blkDimx, size_t blkDimy, size_t grdDimx, size_t grdDimy) {
           dim3 dimBlock( blkDimx, blkDimy );
     	   dim3 dimGrid( grdDimx, grdDimy );
           new_qss<<<dimGrid, dimBlock>>>(kf, kb, MM, Y, Yp, Yc,
	   		      h, alpha, q, p, qp, pp,	
			      q_bar, p_bar, alpha_bar, Ypo, nspec,
			      numreac, numcell, dt_global, debugging);
	}
}
