#pragma OPENCL EXTENSION cl_khr_fp64 : enable

void alpha_func(__global double *p,
                double dt_chem,                            
                __global double *alpha,                    
                int nspec, idx, int numcell)               
{                                                          
double rr;                                                 
for(int i=0; i < nspec; i++)                               
{                                                          
rr = 1.0/(p[idx+numcell*i]*dt_chem+1.0e-50);        
alpha[idx+numcell*i] = (180.0*rr*rr*rr+60.0*rr*rr+11.0*rr+1.0)/(360.0*rr*rr*rr+60.0*rr*rr+12.0*rr+1.0);     
}                                                          
}                                                          

void ydot(__global double *Y,                              
           __global double *MM,                           
           __global double *kf,                           
           __global double *kb,                           
           __global double *p,                            
           __global double *q,                             
           idx, int numcell)                               
{            
<insert-source-terms-here>
}                                                          

__kernel                                                   
void new_qss(__global double *kf,                          
             __global double *kb,                          
             __global double *MM,                          
             __global double *Y,                          
             __global double *Yp,                          
             __global double *Yc,                          
             __global double *h,                           
             __global double *alpha,                       
             __global double *q,                           
             __global double *p,                           
             __global double *qp,                          
             __global double *pp,                          
             __global double *q_bar,                       
             __global double *p_bar,                       
             __global double *alpha_bar,                   
             __global double *Ypo,                         
             const  int nspec,                             
             const  int numreac,                           
             const  int numcell,                           
             double dt_global,                             
             __global double *debugging)                   
{                                                          
private double dt_flow = dt_global;                                     
private double t = 0.0;                                                  
private int flag;                                                                   
private double e = 0.001;                                                           
private double epsilon = e/2.0;                                                      
private int step_fail;                                                             
int idx = get_global_id(0);                                                        
private double dt_chem;                                                              
dt_chem = h[idx];                                                                     
private int max_iterations = 10;                                                     
private int count = 0;                                                               
private double x0 = 0;                                                            
private double cond = 0.0;                                                         
private double sigma = 0.0;                                                        
private double dt_chem_temp = 0;

while (t < dt_flow)                                                  
 {                                                                     
  count = 0;                                                             
  flag = 0;                                                                       
  step_fail = 1;                                                       
  ydot(Y, MM, kf, kb, p, q, idx, numcell);                                        
  alpha_func(p, dt_chem, alpha, nspec, idx, numcell);             
  for (int i = 0; i < nspec; i++)                                                   
       {                                                                        
         Yp[idx+numcell*i] = Y[idx+numcell*i] + (dt_chem*(q[idx+numcell*i]-p[idx+numcell*i]*Y[idx+numcell*i]))/(1.0+alpha[idx+numcell*i]*dt_chem*p[idx+numcell*i]);        
         Ypo[idx+numcell*i] = Yp[idx+numcell*i];      
       }                                                                        

  while (step_fail == 1)                                                        
       {                                                                          
         step_fail = 0;  
         sigma = 1.0e-50;   
         cond = 0.0; 
         ydot(Yp, MM, kf, kb, pp, qp, idx, numcell);      
         for (int j = 0; j < nspec; j++)           
              {            
                p_bar[idx+numcell*j] = 0.5*(pp[idx+numcell*j]+p[idx+numcell*j]); 
              }  
         alpha_func(p_bar, dt_chem, alpha_bar, nspec, idx, numcell); 
         for (int r=0; r < nspec; r++) 
              {    
                q_bar[idx+numcell*r] = alpha_bar[idx+numcell*r]*qp[idx+numcell*r]+(1.0-alpha_bar[idx+numcell*r])*q[idx+numcell*r]; 
              }     
        for (int k = 0; k < nspec; k++) 
              { 
  Yc[idx+numcell*k] = Y[idx+numcell*k] + (dt_chem*(q_bar[idx+numcell*k]-p_bar[idx+numcell*k]*Y[idx+numcell*k]))/(1.0+alpha_bar[idx+numcell*k]*dt_chem*p_bar[idx+numcell*k]);
              }    
         for (int l = 0; l < nspec; l++) 
              { 
                cond = (fabs(Yc[idx+numcell*l]-Ypo[idx+numcell*l])); 
                if (cond >= (e*(Yc[idx+numcell*l]+1.0e-10)))
                    {        
                      step_fail = 1;      
                     }       
              }     
         count = count + 1;    
         if (count > max_iterations) 
             {   
               step_fail = 0; 
               flag  = 1; 
              } 
        if (step_fail == 1)     
             {         
               for (int u = 0; u < nspec; u++)  
                    {   
                      Yp[idx+numcell*u] = Yc[idx+numcell*u];       
                    }      
             }   
       }
 
  if (flag == 0) 
       {   
         for (int m = 0; m < nspec; m++) 
             {
              Y[idx+numcell*m] = Yc[idx+numcell*m];
             }  
          t = t + dt_chem; 
          h[idx] = dt_chem_temp; 
        } 
	
  for (int n=0; n < nspec; n++)
       { 
       cond = (fabs(Yc[idx+numcell*n]-Ypo[idx+numcell*n])/(epsilon*Yc[idx+numcell*n])); 
       if (cond > sigma) 
           {    
             sigma = cond;
           } 
       }

  if (sigma <= 0.0)
     {
      dt_chem_temp = dt_chem;
     }
  else
     {
      x0 = sigma;
      x0 = x0 - 0.5*(x0*x0 - sigma)/x0;
      x0 = x0 - 0.5*(x0*x0 - sigma)/x0;
      x0 = x0 - 0.5*(x0*x0 - sigma)/x0;
      dt_chem_temp = dt_chem * ((1.0/x0) + 0.005);
     }

  if (flag == 0)  
      {   
       if(dt_chem_temp/dt_chem > 1.15)
          {     
            dt_chem = 1.15*dt_chem;   
          }    
       else{     
            dt_chem = dt_chem_temp;        
            }       
      }       

  if (flag == 1)    
      {    
       if(dt_chem_temp/dt_chem > (1.0/3.0))   
          {        
            dt_chem = (1.0/3.0)*dt_chem;    
          }     
       else if(dt_chem_temp/dt_chem < 0.01)
          {    
            dt_chem = 0.01*dt_chem;
           }  
       else
	  {
	   dt_chem = dt_chem_temp;
          }
      }
  
  dt_chem_temp = dt_chem;  
  dt_chem = fmin(dt_flow - t, dt_chem); 
  }    
} 
