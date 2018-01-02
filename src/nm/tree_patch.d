
/**
 * tree_patch.d
 *
 * Contains classes that can be used to build a look-up table (LUT) for a bi-variate function F(x,y). The table can then be searched.
 * The method divides the x,y domain into patches organised in a binomial tree. 
 * Patches are adaptively sized based on the level of error in that patch.
 * The method is based on "Fast Evaluation of Complex Equations of State" by Luke & Collins (2013)
 *
 * Author: Jonathan Ho
 * Date: 17-09-2015
 * Latest Revision: 03-07-2016
 */
module nm.tree_patch;
import std.stdio;
import std.math;
import std.mathspecial;
import std.algorithm; 
import std.string;
import std.conv;
import std.datetime;


alias F_xy = double function(double, double);
alias F_transform = double[2] function(double, double);

static double[16][16] B_inv = 
        [[1.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000],
        [-0.8333333333333335, 3.0000000000000000, -1.5000000000000000, 0.3333333333333333, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000],
        [0.3333333333333334, -1.5000000000000000, 3.0000000000000000, -0.8333333333333333, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000],
        [0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 1.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000],
        [-0.8333333333333335, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 3.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -1.5000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.3333333333333333, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000],
        [0.6944444444444451, -2.5000000000000004, 1.2500000000000002, -0.2777777777777777, -2.5000000000000004, 9.0000000000000000, -4.5000000000000000, 0.9999999999999997, 1.2500000000000002, -4.5000000000000000, 2.2500000000000000, -0.4999999999999997, -0.2777777777777777, 0.9999999999999998, -0.4999999999999997, 0.1111111111111111],
        [-0.2777777777777780, 1.2500000000000002, -2.5000000000000004, 0.6944444444444444, 1.0000000000000002, -4.5000000000000000, 9.0000000000000000, -2.4999999999999996, -0.5000000000000001, 2.2500000000000000, -4.5000000000000000, 1.2499999999999996, 0.1111111111111111, -0.4999999999999998, 0.9999999999999996, -0.2777777777777777],
        [0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.8333333333333335, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 3.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -1.5000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.3333333333333333],
        [0.3333333333333334, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -1.5000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 3.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.8333333333333333, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000],
        [-0.2777777777777780, 1.0000000000000002, -0.5000000000000001, 0.1111111111111111, 1.2500000000000002, -4.5000000000000000, 2.2500000000000000, -0.4999999999999998, -2.5000000000000004, 9.0000000000000000, -4.5000000000000000, 0.9999999999999996, 0.6944444444444444, -2.4999999999999996, 1.2499999999999996, -0.2777777777777778],
        [0.1111111111111112, -0.5000000000000001, 1.0000000000000002, -0.2777777777777778, -0.5000000000000001, 2.2500000000000000, -4.5000000000000000, 1.2499999999999996, 1.0000000000000002, -4.5000000000000000, 9.0000000000000000, -2.4999999999999991, -0.2777777777777778, 1.2499999999999996, -2.4999999999999991, 0.6944444444444443],
        [0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.3333333333333334, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -1.5000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 3.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.8333333333333333],
        [0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 1.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000],
        [0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.8333333333333335, 3.0000000000000000, -1.5000000000000000, 0.3333333333333333],
        [0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.3333333333333334, -1.5000000000000000, 3.0000000000000000, -0.8333333333333333],
        [0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 1.0000000000000000]];
class patch{
        double x_lo;
        double x_hi;
        double y_lo;
        double y_hi;
        double[16] bs;//controlpoints
        
        const double x_m(){return (x_lo + x_hi)*0.5;}
        
        const double y_m(){return (y_lo + y_hi)*0.5;}
                
        double aspectRatio(){return (y_hi-y_lo)/(x_hi-x_lo);}
        
        double area(){return (y_hi-y_lo)*(x_hi-x_lo);}
        
        double transformedArea(F_transform F){
                /*
                p1-------p2
                |         |
                |         |
                |         |
                p0-------p3
                */
                //transforms x, y into some other two-dimensional plane based on F and calculates the area
                //of the quadilateral in this area
                double[2] p0 = F(x_lo,y_lo);
                double[2] p1 = F(x_lo, y_hi);
                double[2] p2 = F(x_hi, y_hi);
                double[2] p3 = F(x_hi, y_lo);
                double[2] A = [p1[0] - p0[0], p1[1] - p0[1]];
                double[2] B = [p2[0] - p0[0], p2[1] - p0[1]];
                double[2] C = [p3[0] - p0[0], p3[1] - p0[1]];
                return 0.5*(fabs(A[0]*B[1]-B[1]*A[0]) + fabs(B[0]*C[1]-B[1]*C[0]));
        }               
        override string toString() const
        {
                return format("This is a patch with x_lo: %s, x_hi %s, y_lo: %s, y_hi: %s", x_lo, x_hi, y_lo, y_hi);
        }
        this(in double x_lo, in double x_hi, in double y_lo, in double y_hi)
        {
                this.x_lo = x_lo;
                this.x_hi = x_hi;
                this.y_lo = y_lo;
                this.y_hi = y_hi;               
        }
        const double interpolateF(double x, double y){
                assert(((x>=x_lo)&&(x<=x_hi)),"x not bounded properly");
                assert(((y>=y_lo)&&(y<=y_hi)),format("y not bounded properly, y: %s is not in the bracket [%s, %s]",y,y_lo,y_hi));
                double u = (x - x_lo)/(x_hi-x_lo);
                double v = (y - y_lo)/(y_hi-y_lo);
                double[4][4][4] bs_k;//now with the k-th dimension
                int i_bs = 0;
                foreach(ref row;bs_k[3]){
                        foreach(ref element; row){
                                element = bs[i_bs];
                                i_bs++;
                                }
                }
                for(int k = 2; k != -1; k--){
                        for(int i = 0; i != k+1; i++){
                                for(int j = 0; j != k+1; j++){
                                        bs_k[k][j][i] = (1-u)*(1-v)*bs_k[k+1][j][i]+u*(1-v)*bs_k[k+1][j][i+1]+(1-u)*v*bs_k[k+1][j+1][i] + u*v*bs_k[k+1][j+1][i+1];
                                        assert(!isNaN(bs_k[k][j][i]),format("i: %s, j: %s, k: %s",i,j,k));
                                }
                        }
                }
                return bs_k[0][0][0];
                }
        
        void getControlPoints(F_xy F)
                {
                double u;
                double v;
                double x;
                double y;
                double[16] fs;
                for(int n = 0; n!=4; n++ ){
                        for(int m = 0; m!=4; m++){//parallelisation possible here
                                u = m/3.0; v = n/3.0;
                                x = (1.0 - u)*this.x_lo + u*this.x_hi;
                                y = (1.0 - v)*this.y_lo + v*this.y_hi;
                                fs[4*n+m] = F(x,y);
                                if (isNaN(fs[4*n+m])) throw new Exception(format("Function calculated NaN at x: %s, y: %s", x,y));
                        }
                }
                foreach(i, ref b; this.bs){
                        b = 0.0;
                        foreach(j, f; fs) b += B_inv[i][j]*f;
                        }
                
                }
        double evalError(F_xy F, double npoints = 10)
        {       
                //EVALUATES ERROR ON A MEAN-SQUARED ERROR METHOD
                //ON A GRID OF npoints x npoints
                double error2 = 0;
                double x;
                double y;
                
                for(int i = 0; i != npoints; i++){
                        for(int j = 0; j != npoints; j++){
                                x = i/npoints*(this.x_hi - this.x_lo) + x_lo;
                                y = j/npoints*(this.y_hi - this.y_lo) + y_lo;
                                error2 += pow(this.interpolateF(x,y)/F(x,y) - 1,2);
                                
                                
                        }
                }
                return sqrt(error2/npoints/npoints);
        }
        double evalMaxError(F_xy F, double npoints = 20)
        {       
                //EVALUATES MAX ERROR BY SAMPLING npoints x npoints
                double maxError = 0;
                double x;
                double y;
                
                for(int i = 0; i != npoints; i++){
                        for(int j = 0; j != npoints; j++){
                                x = i/npoints*(this.x_hi - this.x_lo) + x_lo;
                                y = j/npoints*(this.y_hi - this.y_lo) + y_lo;
                                maxError = max(fabs(this.interpolateF(x,y)/F(x,y) - 1),maxError);
                        }
                }
                return maxError;
        }       
        
        patch splitX_L(){
                double x_m = (x_lo + x_hi)/2.0;
                return new patch(x_lo, x_m, y_lo, y_hi);
        }
        patch splitX_R(){
                double x_m = (x_lo + x_hi)/2.0;
                return new patch(x_m, x_hi, y_lo, y_hi);
                }
        patch splitY_L(){
                double y_m = (y_lo + y_hi)/2.0;
                return new patch(x_lo,x_hi,y_lo,y_m);
        }
        patch splitY_R(){
                double y_m = (y_lo + y_hi)/2.0;
                return new patch(x_lo,x_hi,y_m,y_hi);
        }
        double[5] justOutsideCoordinates(char boundary, double delta){
                switch (boundary){
                        case 'N':
                                return [y_hi + 0.5*delta, 
                                7.0/8.0*x_lo + 1.0/8.0*x_hi, 
                                5.0/8.0*x_lo + 3.0/8.0*x_hi, 
                                3.0/8.0*x_lo + 5.0/8.0*x_hi, 
                                1.0/8.0*x_lo + 7.0/8.0*x_hi];
                        case 'E':
                                return [x_hi + 0.5*delta, 
                                7.0/8.0*y_lo + 1.0/8.0*y_hi, 
                                5.0/8.0*y_lo + 3.0/8.0*y_hi, 
                                3.0/8.0*y_lo + 5.0/8.0*y_hi, 
                                1.0/8.0*y_lo + 7.0/8.0*y_hi];
                        case 'S':
                                return [y_lo - 0.5*delta, 
                                7.0/8.0*x_lo + 1.0/8.0*x_hi, 
                                5.0/8.0*x_lo + 3.0/8.0*x_hi, 
                                3.0/8.0*x_lo + 5.0/8.0*x_hi, 
                                1.0/8.0*x_lo + 7.0/8.0*x_hi];
                        case 'W':
                                return [x_lo - 0.5*delta, 
                                7.0/8.0*y_lo + 1.0/8.0*y_hi, 
                                5.0/8.0*y_lo + 3.0/8.0*y_hi, 
                                3.0/8.0*y_lo + 5.0/8.0*y_hi, 
                                1.0/8.0*y_lo + 7.0/8.0*y_hi];
                        default:
                                throw new Exception("invalid boundary specified for justOutsideCoordinates");
                }
        }
        double[4] justOutsideCorner(int corner, double delta){
                switch (corner){
                        case 0:
                                return [x_lo - 0.5*delta,
                                x_lo + 0.5*delta, 
                                y_lo - 0.5*delta, 
                                y_lo + 0.5*delta];
                        case 1:
                                return [x_hi + 0.5*delta,
                                x_hi - 0.5*delta, 
                                y_lo - 0.5*delta, 
                                y_lo + 0.5*delta];
                        case 2:
                                return [x_hi + 0.5*delta,
                                x_hi - 0.5*delta, 
                                y_hi + 0.5*delta, 
                                y_hi - 0.5*delta];
                        case 3:
                                return [x_lo - 0.5*delta,
                                x_lo + 0.5*delta, 
                                y_hi + 0.5*delta, 
                                y_hi - 0.5*delta];
                        default:
                                throw new Exception("invalid boundary specified for justOutsideCoordinates");
                }
        }
        const double[3] cornerDerivatives(int vertex){
                //return f_x , f_y, f_xy at patch corners
                //Denote vertices as
                //  3----2
                //  |    |
                //  0----1
                double f_x;
                double f_y; 
                double f_xy;
                switch (vertex){
                        case 0:
                                f_x = 3.0*(bs[1] - bs[0])/(x_hi - x_lo);
                                f_y = 3.0*(bs[4] - bs[0])/(y_hi - y_lo);
                                f_xy = 9.0*(bs[5]+bs[0] - bs[1] - bs[4])/(x_hi - x_lo)/(y_hi - y_lo);
                                return [f_x, f_y, f_xy];
                        case 1:
                                f_x = 3.0*(bs[3] - bs[2])/(x_hi - x_lo);
                                f_y = 3.0*(bs[7] - bs[3])/(y_hi - y_lo);
                                f_xy = 9.0*(-bs[3] - bs[6] + bs[2] + bs[7])/(x_hi - x_lo)/(y_hi - y_lo);
                                return [f_x, f_y, f_xy];
                        case 2:
                                f_x = 3.0*(bs[15] - bs[14])/(x_hi - x_lo);
                                f_y = 3.0*(bs[15] - bs[11])/(y_hi - y_lo);
                                f_xy = 9.0*(bs[15]+bs[10] - bs[11] - bs[14])/(x_hi - x_lo)/(y_hi - y_lo);
                                return [f_x, f_y, f_xy];
                        case 3:
                                f_x = 3.0*(bs[13] - bs[12])/(x_hi - x_lo);
                                f_y = 3.0*(bs[12] - bs[8])/(y_hi - y_lo);
                                f_xy = 9.0*(-bs[12] - bs[9] + bs[8] + bs[13])/(x_hi - x_lo)/(y_hi - y_lo);
                                return [f_x, f_y, f_xy];
                        default:
                                throw new Exception("invalid vertex requested for corner derivatives");
                }
        }
        const double[3] midsideDerivatives(char side){
                double[4] alphas = [-1.0, -3.0, 3.0, 1.0];
                double[4] betas = [1.0, 3.0, 3.0, 1.0];
                int[4] idxs = [0, 1, 2, 3];
                double f_x = 0.0; double f_y = 0.0; double f_xy = 0.0;
                int x_sign;
                int y_sign;
                double[16] rotatedBs = rotateControlPoints(side);
                switch (side){
                        case 'N':
                                x_sign = -1; y_sign = -1;
                                break;
                        case 'E':
                                x_sign = 1; y_sign = -1;
                                break;
                        case 'S': 
                                x_sign = 1; y_sign = 1;
                                break;
                        case 'W':
                                x_sign = -1; y_sign = 1;
                                break;
                        default:
                                throw new Exception("Invalid char specified for side");
                        }
                foreach(i; idxs){
                                f_x += alphas[i]*rotatedBs[i];
                                f_y += betas[i]*(rotatedBs[i+4] - rotatedBs[i]);
                                f_xy += alphas[i]*(rotatedBs[i +4] - rotatedBs[i]);
                        }
                f_x = x_sign*f_x*0.5;
                f_y = y_sign*f_y*0.5*0.5*0.5*3;
                f_xy = x_sign*y_sign*f_xy*0.5*3/(x_hi - x_lo)/(y_hi - y_lo);
                if (side=='N'||side=='S') return [f_x/(x_hi - x_lo), f_y/(y_hi - y_lo), f_xy];
                        else return [f_y/(x_hi - x_lo), f_x/(y_hi - y_lo), f_xy];
        }
        const double cornerValue(int vertex){
                //return f at patch corners
                //Denote vertices as
                //  3----2
                //  |    |
                //  0----1
                switch (vertex){
                        case 0:
                                return this.bs[0];
                        case 1:
                                return this.bs[3];
                        case 2: 
                                return this.bs[15];
                        case 3:
                                return this.bs[12];
                        default:
                                throw new Exception("invalid vertex specified");
                }
        }
        const double[16] rotateControlPoints(char side){
                //reorders control points so that the patch appears rotated
                double[16] rotatedBs;
                switch (side){
                        case 'N':
                                for(int i = 0; i != 16; i++){
                                        rotatedBs[i] = this.bs[15-i];
                                }
                                return rotatedBs;
                        case 'E':
                                int[16] rotateidx = [3, 7, 11, 15, 2, 6, 10, 14, 1, 5, 9, 13, 0, 4, 8, 12]; //this was PDiddy's Idea
                                for(int i = 0; i != 16; i++){
                                        rotatedBs[i] = this.bs[rotateidx[i]];
                                }
                                return rotatedBs;
                        case 'S':
                                return this.bs;
                        case 'W':
                                int[16] rotateidx = [12, 8, 4, 0, 13, 9, 5, 1, 14, 10, 6, 2, 15, 11, 7, 3];
                                for(int i = 0; i != 16; i++){
                                        rotatedBs[i] = this.bs[rotateidx[i]];
                                }
                                return rotatedBs;
                        default:
                                throw new Exception("Invalid char specified for side");
                }
        }
        double rewriteControlPoints(int vertex, double f, double f_x, double f_y, double f_xy){
                //rewrite control points based on value of function at vertex and derivatives
                //  3----2
                //  |    |
                //  0----1
                auto oldbs = bs;
                switch (vertex){
                        case 0:
                                bs[0] = f;//0,0
                                bs[1] = f_x*(x_hi - x_lo)/3.0 + bs[0];//0,1
                                bs[4] = f_y*(y_hi - y_lo)/3.0 + bs[0];//1,0
                                bs[5] = f_xy*(x_hi - x_lo)*(y_hi - y_lo)/9.0 - bs[0] + bs[1] + bs[4];//1,1
                                break;
                        case 1:
                                bs[3] = f;
                                bs[2] = -f_x*(x_hi - x_lo)/3.0 + bs[3];
                                bs[7] = f_y*(y_hi - y_lo)/3.0 + bs[3];
                                bs[6] = -f_xy*(x_hi - x_lo)*(y_hi - y_lo)/9.0 - bs[3] + bs[2] + bs[7]; 
                                break;
                        case 2:
                                bs[15] = f;
                                bs[14] = -f_x*(x_hi - x_lo)/3.0 + bs[15];
                                bs[11] = -f_y*(y_hi - y_lo)/3.0 + bs[15];
                                bs[10] = f_xy*(x_hi - x_lo)*(y_hi - y_lo)/9.0 - bs[15] + bs[14] + bs[11];
                                break;
                        case 3:
                                bs[12] = f;
                                bs[13] = f_x*(x_hi - x_lo)/3.0 + bs[12];
                                bs[8] = -f_y*(y_hi - y_lo)/3.0 + bs[12];
                                bs[9] = -f_xy*(x_hi - x_lo)*(y_hi - y_lo)/9.0 - bs[12] + bs[8] + bs[13];
                                break;
                        //case 3:
                        //      {};
                        default:
                                throw new Exception("invalid vertex");
                }
                double maxDiff = 0;
                //foreach (i, ref diff; diffBs) diff = fabs(oldbs[i]/bs[i]-1.0);
                for(int i = 0; i != 16; i++) maxDiff = max(maxDiff, fabs(oldbs[i]/bs[i] - 1.0));
                return maxDiff;
        }
} 


class TreeNode {
        patch nodePatch;
        TreeNode* left;
        TreeNode* right;
        TreeNode* parent;
        int idx;
        char splitID;
        int[] bigNeighbours;
        int[] smallNeighbours;
        
        this(patch myPatch){
                this.nodePatch = myPatch;}
        override string toString() const 
        {
                return format("Node with idx: %s, splitID: %s, left_idx: %s, right_idx: %s", idx, splitID, left, right);
        }
        void writeData(string filename = "Tree.dat",char lastFlag = 'N')
        {
                File tFile = File(filename, "a");
                tFile.writeln(idx);
                if (splitID != 'N') {
                        tFile.writeln((*left).idx);
                        tFile.writeln((*right).idx);
                        }
                else{
                        tFile.writeln(0);
                        tFile.writeln(0);
                        }
                tFile.writeln(splitID);
                tFile.writefln("%.16f",nodePatch.x_lo);
                tFile.writefln("%.16f",nodePatch.x_hi);
                tFile.writefln("%.16f",nodePatch.y_lo);
                tFile.writefln("%.16f",nodePatch.y_hi);
                foreach(b; nodePatch.bs) tFile.writef("%.16e ", b);
                if (lastFlag == 'N') tFile.write("\n");
                                
        }
}
class Tree {
        TreeNode[] Nodes;//the actual tree
        double x_lo;//overall patch that the Tree covers
        double x_hi;
        double y_lo;
        double y_hi;
        //OPTIONAL FOR STORING THE ORIGINAL BOUNDS THAT THE TREE WAS TRANSFORMED FROM
        double X_min;
        double X_max;
        double Y_min;
        double Y_max;
        double globalAspectRatio;
        double globalMinArea;
        double globalMaxError;
        int refinedFlag = 0;
        
        //----------------------------------------------------------
        //constructor which fills the tree with its first node
        this(double x_lo, double x_hi, double y_lo, double y_hi)
        {
                this.globalAspectRatio = (y_hi-y_lo)/(x_hi - x_lo);
                this.globalMinArea = 1e7; 
                this.globalMaxError = 0.001;
                this.x_lo = x_lo;
                this.x_hi = x_hi;
                this.y_lo = y_lo;
                this.y_hi = y_hi;
                Nodes ~= new TreeNode(new patch(x_lo, x_hi, y_lo, y_hi));//append the first patch
                
                }
        
        void grow(F_xy F, F_transform F_t, TreeNode CurrentTreeNode){
        //grows the tree based on the last node
                int idx = to!int(Nodes.length-1);
                //writefln("-------------------------idx: %s -----------------------------", idx);
                CurrentTreeNode.idx = idx;
                if (Nodes.length == 0){throw new Exception("There are no nodes in the tree");}
                CurrentTreeNode.nodePatch.getControlPoints(F);//get control points to evaluate error
                if((CurrentTreeNode.nodePatch.evalError(F) > globalMaxError)&(CurrentTreeNode.nodePatch.transformedArea(F_t) > globalMinArea)){
                        if (CurrentTreeNode.nodePatch.aspectRatio() < globalAspectRatio){
                                CurrentTreeNode.splitID = 'X';
                                this.Nodes ~= new TreeNode(CurrentTreeNode.nodePatch.splitX_L);
                                CurrentTreeNode.left = &Nodes[$-1]; 
                                this.grow(F,F_t,Nodes[$-1]);
                                
                                this.Nodes ~= new TreeNode(CurrentTreeNode.nodePatch.splitX_R);
                                CurrentTreeNode.right = &Nodes[$-1];
                                this.grow(F,F_t,Nodes[$-1]);
                        }
                        else{
                                CurrentTreeNode.splitID = 'Y';
                                this.Nodes ~= new TreeNode(CurrentTreeNode.nodePatch.splitY_L);
                                CurrentTreeNode.left = &Nodes[$-1];
                                this.grow(F,F_t,Nodes[$-1]);
                                this.Nodes ~= new TreeNode(CurrentTreeNode.nodePatch.splitY_R);
                                CurrentTreeNode.right = &Nodes[$-1];
                                this.grow(F,F_t,Nodes[$-1]);
                                }
                        }
                else{
                        CurrentTreeNode.splitID = 'N';
                        }
        }
        const const(patch) search(double x, double y){
        //returns the patch to the node in the tree that bounds x and y
        //had to use pointers as it was const
        const(TreeNode)* currentNode = &Nodes[0];//start at first node
                if ((x<currentNode.nodePatch.x_lo)||(x>currentNode.nodePatch.x_hi)){
                        throw new Exception(format("x (probably re-parameterized internal energy) not bounded properly by look up table, x: %s is not in the bracket [%s, %s]",x,currentNode.nodePatch.x_lo,currentNode.nodePatch.x_hi));
                }
                if ((y<currentNode.nodePatch.y_lo)||(y>currentNode.nodePatch.y_hi)){
                        throw new Exception(format("y (probably re-parametrized density) not bounded properly by look up table, y: %s is not in the bracket [%s, %s]",y,currentNode.nodePatch.y_lo,currentNode.nodePatch.y_hi));
                }
                while((*currentNode).splitID != 'N') {
                        if (currentNode.splitID == 'X'){
                                if (x < (*currentNode).nodePatch.x_m()) currentNode = (*currentNode).left;
                                else currentNode = (*currentNode).right;
                                }
                        else if (currentNode.splitID == 'Y'){
                                if (y < (*currentNode).nodePatch.y_m()) currentNode = (*currentNode).left;
                                else currentNode = (*currentNode).right;
                                }
                        
                        } 
                return (*currentNode).nodePatch;
        }
        int searchForNodeID(double x, double y){
        //returns index of the node of the patch in the tree that bounds x and y
        //had to use pointers as it was const
        const(TreeNode)* currentNode = &Nodes[0];//start at first node
                if ((x<currentNode.nodePatch.x_lo)||(x>currentNode.nodePatch.x_hi)){
                        throw new Exception(format("x (probably re-parameterized internal energy) not bounded properly by look up table, x: %s is not in the bracket [%s, %s]",x,currentNode.nodePatch.x_lo,currentNode.nodePatch.x_hi));
                }
                if ((y<currentNode.nodePatch.y_lo)||(y>currentNode.nodePatch.y_hi)){
                        throw new Exception(format("y (probably re-parametrized density) not bounded properly by look up table, y: %s is not in the bracket [%s, %s]",y,currentNode.nodePatch.y_lo,currentNode.nodePatch.y_hi));
                }
                while((*currentNode).splitID != 'N') {
                        if (currentNode.splitID == 'X'){
                                if (x < (*currentNode).nodePatch.x_m()) currentNode = (*currentNode).left;
                                else currentNode = (*currentNode).right;
                                }
                        else if (currentNode.splitID == 'Y'){
                                if (y < (*currentNode).nodePatch.y_m()) currentNode = (*currentNode).left;
                                else currentNode = (*currentNode).right;
                                }
                        
                        } 
                return (*currentNode).idx;
        }
        
        void writeLeaves(){
                //write's co-ordinates of all leaf node patches for plotting in python
                File patchFile = File("patch_xy.dat", "w");
                foreach (node; Nodes){
                        if (node.splitID == 'N'){
                                patchFile.writeln([node.nodePatch.x_lo, node.nodePatch.x_hi,node.nodePatch.y_lo,node.nodePatch.y_hi]);
                                }
                        }
                
                }
        void writeLUT(string filename = "Tree.dat"){
                File tFile = File(filename, "w");
                tFile.writeln(globalMaxError);
                tFile.writeln(globalMinArea);
                tFile.writeln(X_min);
                tFile.writeln(X_max);
                tFile.writeln(Y_min);
                tFile.writeln(Y_max);
                tFile.close();
                foreach(i, Node; Nodes) {
                        if(i != Nodes.length - 1) Node.writeData(filename);
                        else Node.writeData(filename,'L');
                }
                }
        double minDeltaX(){
                double minDeltaX = X_max - X_min;
                foreach(node; Nodes){
                        minDeltaX = min(minDeltaX,node.nodePatch.x_hi - node.nodePatch.x_lo);
                }
                return minDeltaX;
        }
        double minDeltaY(){
                double minDeltaY = Y_max - Y_min;
                foreach(node; Nodes){
                        minDeltaY = min(minDeltaY,node.nodePatch.y_hi - node.nodePatch.y_lo);
                }
                return minDeltaY;
        }
        int refine(F_xy F, F_transform F_t){
                //returns the number of patches split for refinement purposes
                double minDeltaXorY = min(this.minDeltaX(),this.minDeltaY());
                double[5] coords;
                int[4] nodeIds;
                int numberRefined = 0;
                foreach (node; Nodes){
                        if (node.splitID == 'N'){
                                foreach(boundary; ['N', 'E', 'S', 'W']){
                                        coords = node.nodePatch.justOutsideCoordinates(boundary,minDeltaXorY);
                                        try{
                                                foreach(i;[0,1,2,3]){
                                                        if ((boundary=='N')|(boundary=='S')) nodeIds[i] = this.searchForNodeID(coords[i+1],coords[0]);
                                                                else nodeIds[i] = this.searchForNodeID(coords[0],coords[i+1]);
                                                }
                                                if ((nodeIds[0]!=nodeIds[1])|(nodeIds[2]!=nodeIds[3])){
                                                        numberRefined += 1;
                                                        if ((boundary=='N')|(boundary=='S')){
                                                                //do splitting and break
                                                                node.splitID = 'X';
                                                                this.Nodes ~= new TreeNode(node.nodePatch.splitX_L);
                                                                node.left = &Nodes[$-1];
                                                                int idx = to!int(Nodes.length-1);
                                                                Nodes[$-1].idx = idx;
                                                                Nodes[$-1].nodePatch.getControlPoints(F);
                                                                Nodes[$-1].splitID = 'N';
                                                                //this.grow(F,F_t,Nodes[$-1]);
                                                                this.Nodes ~= new TreeNode(node.nodePatch.splitX_R);
                                                                node.right = &Nodes[$-1];
                                                                //this.grow(F,F_t,Nodes[$-1]);
                                                                idx = to!int(Nodes.length-1);
                                                                Nodes[$-1].idx = idx;
                                                                Nodes[$-1].nodePatch.getControlPoints(F);
                                                                Nodes[$-1].splitID = 'N';

                                                        }
                                                        else{
                                                                node.splitID = 'Y';
                                                                this.Nodes ~= new TreeNode(node.nodePatch.splitY_L);
                                                                node.left = &Nodes[$-1];
                                                                int idx = to!int(Nodes.length-1);
                                                                Nodes[$-1].idx = idx;
                                                                Nodes[$-1].nodePatch.getControlPoints(F);
                                                                Nodes[$-1].splitID = 'N';
                                                                //this.grow(F,F_t,Nodes[$-1]);
                                                                this.Nodes ~= new TreeNode(node.nodePatch.splitY_R);
                                                                node.right = &Nodes[$-1];
                                                                //this.grow(F,F_t,Nodes[$-1]);
                                                                idx = to!int(Nodes.length-1);
                                                                Nodes[$-1].idx = idx;
                                                                Nodes[$-1].nodePatch.getControlPoints(F);
                                                                Nodes[$-1].splitID = 'N';
                                                        }
                                                        break;
                                                }//end if nodes are different
                                        }//end try
                                        catch (Exception e) {
                                            writeln(e.msg);
                                            //writeln("Tried to refine outside the table, continue on");
                                        }
                                }//end foreach boundary
                        }
                }//end foreach node loop
                return numberRefined;
        }//end refine function
        void recordDependencies(){
                //returns the number of patches split for refinement purposes
                if (this.refinedFlag != 1) throw new Exception("Tree is not refined yet, cannot record dependencies");
                double minDeltaXorY = min(this.minDeltaX(),this.minDeltaY());
                double[5] coords;
                int nodeID1;
                int nodeID2;
                foreach (node; Nodes){
                        if (node.splitID == 'N'){
                                foreach(boundary; ['N', 'E', 'S', 'W']){
                                        coords = node.nodePatch.justOutsideCoordinates(boundary,minDeltaXorY);
                                        try{
                                                        if ((boundary=='N')|(boundary=='S')) {
                                                                nodeID1 = this.searchForNodeID(coords[1],coords[0]);
                                                                nodeID2 = this.searchForNodeID(coords[3],coords[0]);
                                                        }                                                               
                                                        else{
                                                                nodeID1 = this.searchForNodeID(coords[0],coords[1]);
                                                                nodeID2 = this.searchForNodeID(coords[0],coords[3]);
                                                                } 
                                                        if (nodeID1 != nodeID2){//two different patchees
                                                                node.smallNeighbours ~= [nodeID1, nodeID2];
                                                                this.Nodes[nodeID1].bigNeighbours ~= node.idx;
                                                                this.Nodes[nodeID2].bigNeighbours ~= node.idx;
                                                        }

                                                }//end try
                                        catch (Exception e) {
                                            writeln(e.msg);
                                            //writeln("Tried to refine outside the table, continue on");
                                        }
                                }//end foreach boundary
                        }
                }//end foreach node loop
        }//end recordDependencies function
        int[] dependencyOrder(){
                //records order of nodes in which to update
                int[] nodeIDorder;
                int numberAdded=1;
                while (numberAdded){
                        numberAdded = 0;
                        foreach(node; Nodes){
                                if (node.splitID == 'N'){
                                        if (!node.bigNeighbours&&!nodeIDorder.canFind(node.idx)) {
                                                nodeIDorder ~= node.idx;
                                                numberAdded += 1;
                                                foreach(smallNeighbour;node.smallNeighbours) {//remove the bigNeighbour from list
                                                        int[] newBigNeighbours;
                                                        foreach(bigNeighbour; Nodes[smallNeighbour].bigNeighbours) if (bigNeighbour!=node.idx) newBigNeighbours ~= bigNeighbour;
                                                        Nodes[smallNeighbour].bigNeighbours=newBigNeighbours;
                                                }
                                        }//end if the node has no big neighbours
                                }//end if node is leaf
                        }//end foreach node
                }//end while
                return nodeIDorder;
        }
        void makeMidsideContinuous(int[] dependencyOrder){
                //Goes through each leaf node in the tree and checks for an adjacent side which has two patches
                //It then evaluates the function and three derivatives f_x, f_y, f_xy to propagate to the corners of these
                //two patches
                //Checks through leaf nodes in the order specified by dependency Order
                //The following class methods should be used before makeMidsideContinuous:
                //1. refine: splits patches to ensure that each patch has no more than two neighbouring patches on one side
                //2. recordDependencies: records the ids of neighbouring patches depending on whether they are bigger or smaller than the patch
                //3. dependencyOrder: records the dependency of the patches so they are updated in the right order
                //4. run makeMidsideContinuous

                double f;
                double[3] derivs;
                int nodeID1;
                int nodeID2;
                double[5] coords;
                double minDeltaXorY = min(this.minDeltaX(),this.minDeltaY());
                foreach(nodeID;dependencyOrder){
                        double x_lo = Nodes[nodeID].nodePatch.x_lo;
                        double x_hi = Nodes[nodeID].nodePatch.x_hi;
                        double y_lo = Nodes[nodeID].nodePatch.y_lo;
                        double y_hi = Nodes[nodeID].nodePatch.y_hi;

                        if (Nodes[nodeID].splitID == 'N'){
                                foreach(boundary; ['N', 'E', 'S', 'W']){
                                        coords = Nodes[nodeID].nodePatch.justOutsideCoordinates(boundary,minDeltaXorY);
                                        try{
                                                if ((boundary=='N')|(boundary=='S')) {
                                                        nodeID1 = this.searchForNodeID(coords[1],coords[0]);
                                                        nodeID2 = this.searchForNodeID(coords[3],coords[0]);
                                                }                                                               
                                                else{
                                                        nodeID1 = this.searchForNodeID(coords[0],coords[1]);
                                                        nodeID2 = this.searchForNodeID(coords[0],coords[3]);
                                                        } 
                                                }//end try
                                        catch (Exception e) {
                                            writeln(e.msg);
                                            nodeID1=1;nodeID2=1;//writeln("Tried to refine outside the table, continue on");
                                        }
                                                if (nodeID1 != nodeID2){//two different patchees
                                                        double percentagediff;
                                                        double maxDiff=0.01;
                                                        switch(boundary){
                                                                case 'N':
                                                                        //top  left corner
                                                                        derivs = Nodes[nodeID].nodePatch.cornerDerivatives(3);
                                                                        f = Nodes[nodeID].nodePatch.cornerValue(3);
                                                                        Nodes[nodeID1].nodePatch.rewriteControlPoints(0,f,derivs[0],derivs[1],derivs[2]);
                                                                        //MIDDLE TOP
                                                                        derivs = Nodes[nodeID].nodePatch.midsideDerivatives('N');
                                                                        f = Nodes[nodeID].nodePatch.interpolateF(0.5*(x_lo+x_hi), y_hi);
                                                                        Nodes[nodeID1].nodePatch.rewriteControlPoints(1,f,derivs[0],derivs[1],derivs[2]);
                                                                        Nodes[nodeID2].nodePatch.rewriteControlPoints(0,f,derivs[0],derivs[1],derivs[2]);
                                                                        //top right
                                                                        derivs = Nodes[nodeID].nodePatch.cornerDerivatives(2);
                                                                        f = Nodes[nodeID].nodePatch.cornerValue(2);
                                                                        Nodes[nodeID2].nodePatch.rewriteControlPoints(1,f,derivs[0],derivs[1],derivs[2]);
                                                                        break;
                                                                case 'E':
                                                                        //bottom  right corner
                                                                        derivs = Nodes[nodeID].nodePatch.cornerDerivatives(1);
                                                                        f = Nodes[nodeID].nodePatch.cornerValue(1);
                                                                        Nodes[nodeID1].nodePatch.rewriteControlPoints(0,f,derivs[0],derivs[1],derivs[2]);
                                                                        //MIDDLE left
                                                                        derivs = Nodes[nodeID].nodePatch.midsideDerivatives('E');
                                                                        f = Nodes[nodeID].nodePatch.interpolateF(x_hi, 0.5*(y_lo+y_hi));
                                                                        Nodes[nodeID1].nodePatch.rewriteControlPoints(3,f,derivs[0],derivs[1],derivs[2]);
                                                                        Nodes[nodeID2].nodePatch.rewriteControlPoints(0,f,derivs[0],derivs[1],derivs[2]);
                                                                        //top right
                                                                        derivs = Nodes[nodeID].nodePatch.cornerDerivatives(2);
                                                                        f = Nodes[nodeID].nodePatch.cornerValue(2);
                                                                        Nodes[nodeID2].nodePatch.rewriteControlPoints(3,f,derivs[0],derivs[1],derivs[2]);
                                                                        break;
                                                                case 's':
                                                                        //bottom left corner
                                                                        derivs = Nodes[nodeID].nodePatch.cornerDerivatives(0);
                                                                        f = Nodes[nodeID].nodePatch.cornerValue(0);
                                                                        Nodes[nodeID1].nodePatch.rewriteControlPoints(3,f,derivs[0],derivs[1],derivs[2]);
                                                                        //MIDDLE bottom
                                                                        derivs = Nodes[nodeID].nodePatch.midsideDerivatives('S');
                                                                        f = Nodes[nodeID].nodePatch.interpolateF(0.5*(x_lo+x_hi), y_lo);
                                                                        Nodes[nodeID1].nodePatch.rewriteControlPoints(2,f,derivs[0],derivs[1],derivs[2]);
                                                                        Nodes[nodeID2].nodePatch.rewriteControlPoints(3,f,derivs[0],derivs[1],derivs[2]);
                                                                        //bottom right
                                                                        derivs = Nodes[nodeID].nodePatch.cornerDerivatives(1);
                                                                        f = Nodes[nodeID].nodePatch.cornerValue(1);
                                                                        Nodes[nodeID2].nodePatch.rewriteControlPoints(2,f,derivs[0],derivs[1],derivs[2]);
                                                                        break;
                                                                case 'W':
                                                                        //bottom  left corner
                                                                        derivs = Nodes[nodeID].nodePatch.cornerDerivatives(0);
                                                                        f = Nodes[nodeID].nodePatch.cornerValue(0);
                                                                        Nodes[nodeID1].nodePatch.rewriteControlPoints(1,f,derivs[0],derivs[1],derivs[2]);
                                                                        //MIDDLE left
                                                                        derivs = Nodes[nodeID].nodePatch.midsideDerivatives('W');
                                                                        f = Nodes[nodeID].nodePatch.interpolateF(x_lo, 0.5*(y_lo+y_hi));
                                                                        Nodes[nodeID1].nodePatch.rewriteControlPoints(2,f,derivs[0],derivs[1],derivs[2]);
                                                                        Nodes[nodeID2].nodePatch.rewriteControlPoints(1,f,derivs[0],derivs[1],derivs[2]);
                                                                        //top left
                                                                        derivs = Nodes[nodeID].nodePatch.cornerDerivatives(3);
                                                                        f = Nodes[nodeID].nodePatch.cornerValue(3);
                                                                        Nodes[nodeID2].nodePatch.rewriteControlPoints(2,f,derivs[0],derivs[1],derivs[2]);
                                                                        break;
                                                                default:
                                                                        break;
                                                        }//end switch
                                                }//end if one node equals the other


                                                
                                }//end foreach boundary
                        }
                }//end for loop
        }//end make Coninuous function
        void makeVertexContinuous(){
                double f;
                double[3][4] derivs;
                double[3] derivsmins;
                int nodeID1;
                int nodeID2;
                int nodeID3;
                double[4] coords;
                double minDeltaXorY = min(this.minDeltaX(),this.minDeltaY());
                foreach(node;Nodes){
                        double x_lo = node.nodePatch.x_lo;
                        double x_hi = node.nodePatch.x_hi;
                        double y_lo = node.nodePatch.y_lo;
                        double y_hi = node.nodePatch.y_hi;

                        if (node.splitID == 'N'){
                                foreach(vertex; [0, 1, 2, 3]){
                                        coords = node.nodePatch.justOutsideCorner(vertex,minDeltaXorY);
                                        try{
                                                        nodeID1 = this.searchForNodeID(coords[0],coords[3]);
                                                        nodeID2 = this.searchForNodeID(coords[0],coords[2]);
                                                        nodeID3 = this.searchForNodeID(coords[1],coords[2]);
                                                }//end try
                                        catch (Exception e) {
                                            writeln(e.msg);
                                            nodeID1=1;nodeID2=1;
                                        }
                                                if ((nodeID1 != nodeID2)&&(nodeID2 != nodeID3)&&(nodeID1!= nodeID3)){//corner
                                                        switch(vertex){
                                                                case 0:
                                                                        //top  left corner
                                                                        derivs[0] = node.nodePatch.cornerDerivatives(0);
                                                                        derivs[1] = Nodes[nodeID1].nodePatch.cornerDerivatives(1);//need to check these
                                                                        derivs[2] = Nodes[nodeID2].nodePatch.cornerDerivatives(2);
                                                                        derivs[3] = Nodes[nodeID3].nodePatch.cornerDerivatives(3);
                                                                        f = node.nodePatch.cornerValue(0);//potential to smooth this process out here by just picking the control point
                                                                        derivsmins[0] = min(derivs[0][0], derivs[1][0], derivs[2][0]);
                                                                        derivsmins[1] = min(derivs[0][1], derivs[1][1], derivs[2][1]);
                                                                        derivsmins[2] = min(derivs[0][2], derivs[1][2], derivs[2][2]);
                                                                        node.nodePatch.rewriteControlPoints(0,f,derivsmins[0],derivsmins[1],derivsmins[2]);
                                                                        Nodes[nodeID1].nodePatch.rewriteControlPoints(1,f,derivsmins[0],derivsmins[1],derivsmins[2]);
                                                                        Nodes[nodeID2].nodePatch.rewriteControlPoints(2,f,derivsmins[0],derivsmins[1],derivsmins[2]);
                                                                        Nodes[nodeID3].nodePatch.rewriteControlPoints(3,f,derivsmins[0],derivsmins[1],derivsmins[2]);
                                                                        break;
                                                                case 1:
                                                                        //bottom  right corner
                                                                        derivs[0] = node.nodePatch.cornerDerivatives(1);
                                                                        derivs[1] = Nodes[nodeID1].nodePatch.cornerDerivatives(0);//need to check these
                                                                        derivs[2] = Nodes[nodeID2].nodePatch.cornerDerivatives(3);
                                                                        derivs[3] = Nodes[nodeID3].nodePatch.cornerDerivatives(2);
                                                                        f = node.nodePatch.cornerValue(1);//potential to smooth this process out here by just picking the control point
                                                                        derivsmins[0] = min(derivs[0][0], derivs[1][0], derivs[2][0]);
                                                                        derivsmins[1] = min(derivs[0][1], derivs[1][1], derivs[2][1]);
                                                                        derivsmins[2] = min(derivs[0][2], derivs[1][2], derivs[2][2]);
                                                                        node.nodePatch.rewriteControlPoints(1,f,derivsmins[0],derivsmins[1],derivsmins[2]);
                                                                        Nodes[nodeID1].nodePatch.rewriteControlPoints(0,f,derivsmins[0],derivsmins[1],derivsmins[2]);
                                                                        Nodes[nodeID2].nodePatch.rewriteControlPoints(3,f,derivsmins[0],derivsmins[1],derivsmins[2]);
                                                                        Nodes[nodeID3].nodePatch.rewriteControlPoints(2,f,derivsmins[0],derivsmins[1],derivsmins[2]);
                                                                        break;
                                                                case 2:
                                                                        //upper right corner
                                                                        derivs[0] = node.nodePatch.cornerDerivatives(2);
                                                                        derivs[1] = Nodes[nodeID1].nodePatch.cornerDerivatives(3);//need to check these
                                                                        derivs[2] = Nodes[nodeID2].nodePatch.cornerDerivatives(0);
                                                                        derivs[3] = Nodes[nodeID3].nodePatch.cornerDerivatives(1);
                                                                        f = node.nodePatch.cornerValue(2);//potential to smooth this process out here by just picking the control point
                                                                        derivsmins[0] = min(derivs[0][0], derivs[1][0], derivs[2][0]);
                                                                        derivsmins[1] = min(derivs[0][1], derivs[1][1], derivs[2][1]);
                                                                        derivsmins[2] = min(derivs[0][2], derivs[1][2], derivs[2][2]);
                                                                        node.nodePatch.rewriteControlPoints(2,f,derivsmins[0],derivsmins[1],derivsmins[2]);
                                                                        Nodes[nodeID1].nodePatch.rewriteControlPoints(3,f,derivsmins[0],derivsmins[1],derivsmins[2]);
                                                                        Nodes[nodeID2].nodePatch.rewriteControlPoints(0,f,derivsmins[0],derivsmins[1],derivsmins[2]);
                                                                        Nodes[nodeID3].nodePatch.rewriteControlPoints(1,f,derivsmins[0],derivsmins[1],derivsmins[2]);
                                                                        break;
                                                                case 3:
                                                                        //upper  left corner
                                                                        derivs[0] = node.nodePatch.cornerDerivatives(3);
                                                                        derivs[1] = Nodes[nodeID1].nodePatch.cornerDerivatives(2);//need to check these
                                                                        derivs[2] = Nodes[nodeID2].nodePatch.cornerDerivatives(1);
                                                                        derivs[3] = Nodes[nodeID3].nodePatch.cornerDerivatives(0);
                                                                        f = node.nodePatch.cornerValue(3);//potential to smooth this process out here by just picking the control point
                                                                        derivsmins[0] = min(derivs[0][0], derivs[1][0], derivs[2][0]);
                                                                        derivsmins[1] = min(derivs[0][1], derivs[1][1], derivs[2][1]);
                                                                        derivsmins[2] = min(derivs[0][2], derivs[1][2], derivs[2][2]);
                                                                        node.nodePatch.rewriteControlPoints(3,f,derivsmins[0],derivsmins[1],derivsmins[2]);
                                                                        Nodes[nodeID1].nodePatch.rewriteControlPoints(2,f,derivsmins[0],derivsmins[1],derivsmins[2]);
                                                                        Nodes[nodeID2].nodePatch.rewriteControlPoints(1,f,derivsmins[0],derivsmins[1],derivsmins[2]);
                                                                        Nodes[nodeID3].nodePatch.rewriteControlPoints(0,f,derivsmins[0],derivsmins[1],derivsmins[2]);
                                                                        break;
                                                                default:
                                                                        break;
                                                        }//end switch
                                                }//end if one node equals the other


                                                
                                }//end foreach boundary
                        }
                }//end for loop
        }//end make Coninuous function
}//end class Tree

Tree buildTree_fromFile(string filename = "Tree.dat"){
        File treeFile = File(filename,"r");
        int i = 0;
        string[9] lines;
        int idx;
        int[] lefts;
        int[] rights;
        char splitID;
        double x_lo;
        double x_hi;
        double y_lo;
        double y_hi;
        double[16] bs;
        Tree myTree = new Tree(0,100,100,0);
        //First two lines have some tree data
        myTree.globalMaxError = to!double(chomp(treeFile.readln())); 
        myTree.globalMinArea = to!double(chomp(treeFile.readln()));
        myTree.X_min = to!double(chomp(treeFile.readln()));
        myTree.X_max = to!double(chomp(treeFile.readln()));
        myTree.Y_min = to!double(chomp(treeFile.readln()));
        myTree.Y_max = to!double(chomp(treeFile.readln()));
        while (!treeFile.eof()){
                for (int j = 0; j != 9; j++){
                        lines[j] = chomp(treeFile.readln());
                }
                idx = to!int(lines[0]);
                lefts ~= to!int(lines[1]);
                rights ~= to!int(lines[2]);
                splitID = to!char(lines[3]);
                x_lo = to!double(lines[4]);
                x_hi = to!double(lines[5]);
                y_lo = to!double(lines[6]);
                y_hi = to!double(lines[7]);
                foreach(b_i, b; split(lines[8])) bs[b_i] = to!double(b);
                if (i == 0) {
                        myTree.Nodes[0] = new TreeNode(new patch(x_lo,x_hi,y_lo,y_hi));//there is a node already in there
                        myTree.x_lo = x_lo;
                        myTree.x_hi = x_hi;
                        myTree.y_lo = y_lo;
                        myTree.y_hi = y_hi;
                        }
                        else myTree.Nodes ~= new TreeNode(new patch(x_lo,x_hi,y_lo,y_hi));
                myTree.Nodes[i].idx = idx;
                myTree.Nodes[i].splitID = splitID;
                myTree.Nodes[i].nodePatch.bs = bs;
                i++;
        }
        //back fill the pointers now that the whole tree is constructed
        foreach(node_i,node; myTree.Nodes){
                if (node.splitID != 'N') {
                        node.left = &myTree.Nodes[lefts[node_i]];
                        node.right = &myTree.Nodes[rights[node_i]];
                        }
                }
        writefln("Finished reading in EOS look up table from %s", filename);
        return myTree;
}


double[] linspace(double start, double stop, double n){
        double[] ys;
        for(int i = 0; i != n+1;i++){
                ys ~= (stop-start)*i/n + start;
        }
        return ys;
        }




