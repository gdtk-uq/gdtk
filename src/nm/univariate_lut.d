//A module for doing look up on a univariate function
module nm.univariate_lut;
import std.stdio;
import std.math;
import std.mathspecial;
import std.algorithm; 
import std.string;
import std.conv;
import std.datetime;
//quadratic interpolation from:
//http://homepage.math.uiowa.edu/~atkinson/ftp/ENA_Materials/Overheads/sec_4-1.pdf
alias F_x = double function(double);

class uni_lut{
        double[] fs;
        int n;//number of points used
        double xmin;
        double xmax;
        double dx;

        const double lookup(double x){
                //linear lookup
                int i = to!int(trunc(n*(x - xmin)/(xmax - xmin)));
                double p = (x - i*dx - xmin)/dx;
                i = max(i, 0);
                i = min(i, n -1);
                return p*fs[i+1] + (1 - p)*fs[i];
        }
        const double quadlookup(double x){
                //quadratic look up
                int i = to!int(trunc(n*(x - xmin)/(xmax - xmin)));
                double p = (x - i*dx)/dx;
                if (p<0.5) i -= 1;
                i = max(i,0);//don't go too low
                i = min(i, n - 2);//don't go too high
                double x0 = i*dx + xmin; 
                double x1 = (i+1)*dx + xmin; 
                double x2 = (i+2)*dx + xmin;
                double L0 = (x - x1)*(x - x2)/2.0;
                double L1 = -(x - x0)*(x - x2);
                double L2 = (x - x0)*(x - x1)/2.0;
                return (fs[i]*L0 + fs[i+1]*L1 + fs[i+2]*L2)/dx/dx;
                
        }
        void writeLUT(string filename){
                File tFile = File(filename,"w");
                foreach(f; this.fs) tFile.writef("%.16e ", f);
                tFile.write("\n");
                tFile.writeln(n);
                tFile.writeln(xmin);
                tFile.writeln(xmax);
                tFile.writeln(dx);
                writefln("Table written to %s", filename);
        }
        void writeCode(string filename = "uni_lut_literal.d"){
                //writes one line to a textfile that can be used to initialise as local static struct
                File tFile = File(filename,"w");
                tFile.write("uni_lut table_demo = new uni_lut(");
                tFile.writef("%.5f",xmin);
                tFile.write(", ");
                tFile.writef("%.5f",xmax);
                tFile.write(", ");
                tFile.write(n);
                tFile.write(", [");
                foreach(f; this.fs) tFile.writef("%.16e, ", f);
                tFile.write("]); ");
                writefln("Code written to %s", filename);
        }
        this(string filename){
                //constructs from a file that is written previously by writeLUT
                File tFile = File(filename,"r");
                string[5] lines;
                for (int i = 0; i != 5; i++) lines[i] = chomp(tFile.readln());
                foreach(f ; split(lines[0])) this.fs ~= to!double(f);
                this.n = to!int(lines[1]);
                this.xmin = to!double(lines[2]);
                this.xmax = to!double(lines[3]);
                this.dx = to!double(lines[4]);
                writefln("Succesfully constructed table from: %s", filename);
        }
        this(F_x f, double xmin, double xmax, int n){
                this.xmin = xmin;
                this.xmax = xmax;
                this.n = n;
                this.dx = (xmax - xmin)/n;
                for(int i = 0; i != n + 1; i++) {
                        fs ~= f(i*this.dx + xmin);      
                }
        }
        this(double xmin, double xmax, int n, double[] fs) 
        {
                this.xmin = xmin;
                this.xmax = xmax;
                this.n = n;
                this.dx = (xmax - xmin)/n;
                this.fs = fs;
        }       
        

}

