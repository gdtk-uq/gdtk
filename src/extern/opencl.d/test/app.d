import std.stdio;
import cl;

void main()
{
	writeln("Edit source/app.d to start your project.");
	cl_uint a, b;
	cl_platform_id id;
	clGetPlatformIDs(a, &id, &b);
}

