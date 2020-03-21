import std.stdio, std.algorithm, std.typecons;
import std.string, std.conv, std.array, std.utf;
import std.file : readText;

void main(string[] args)
{
    auto mpi = args[1].readText();

    Nullable!(int)[string] intEntries;
    foreach(m; ints)
    {
        intEntries[m] = Nullable!int.init;
    }

    string[string] typeEntries;
    foreach(m; types)
    {
        typeEntries[m] = null;
    }

    foreach(line; mpi.splitLines())
    {
        line = line.byChar.map!(c => c == '\t' ? ' ' : c).to!string;
        auto r = line.findSplitAfter("#define")[1].stripLeft.findSplit(" ");
        r[0] = r[0].stripRight;
        r[2] = r[2].strip;
        if(r[0].empty) continue;
        auto pI = r[0] in intEntries;
        if(pI) *pI = r[2].to!int;
        else
        {
            auto pT = r[0] in typeEntries;
            if(r[2] == "long long") r[2] = "long";
            else if(r[2] == "long") r[2] = "c_long";
            if(pT) *pT = r[2];
        }
    }

    foreach(key, val; intEntries)
    {
        if(val.isNull)
            stderr.writeln("no entry found for " ~ key);
        else
            writeln("enum " ~ key ~ " = " ~ val.to!string ~ ';');
    }
    foreach(key, val; typeEntries)
    {
        if(val is null)
            stderr.writeln("no entry found for " ~ key);
        else
            writeln("alias " ~ key ~ " = " ~ val ~ ';');
    }
}

immutable ints = [
    "OPAL_BUILD_PLATFORM_COMPILER_FAMILYID",
    "OPAL_BUILD_PLATFORM_COMPILER_VERSION",
    "OPAL_STDC_HEADERS",
    "OPAL_HAVE_ATTRIBUTE_DEPRECATED",
    "OPAL_HAVE_ATTRIBUTE_DEPRECATED_ARGUMENT",
    "OPAL_HAVE_SYS_TIME_H",
    "OPAL_HAVE_SYS_SYNCH_H",
    "OPAL_HAVE_LONG_LONG",
    "OPAL_SIZEOF_BOOL",
    "OPAL_SIZEOF_INT",
    "OPAL_MAX_DATAREP_STRING",
    "OPAL_MAX_ERROR_STRING",
    "OPAL_MAX_INFO_KEY",
    "OPAL_MAX_INFO_VAL",
    "OPAL_MAX_OBJECT_NAME",
    "OPAL_MAX_PORT_NAME",
    "OPAL_MAX_PROCESSOR_NAME",
    "OMPI_HAVE_FORTRAN_LOGICAL1",
    "OMPI_HAVE_FORTRAN_LOGICAL2",
    "OMPI_HAVE_FORTRAN_LOGICAL4",
    "OMPI_HAVE_FORTRAN_LOGICAL8",
    "OMPI_HAVE_FORTRAN_INTEGER1",
    "OMPI_HAVE_FORTRAN_INTEGER16",
    "OMPI_HAVE_FORTRAN_INTEGER2",
    "OMPI_HAVE_FORTRAN_INTEGER4",
    "OMPI_HAVE_FORTRAN_INTEGER8",
    "OMPI_HAVE_FORTRAN_REAL16",
    "OMPI_HAVE_FORTRAN_REAL2",
    "OMPI_HAVE_FORTRAN_REAL4",
    "OMPI_HAVE_FORTRAN_REAL8",
    "HAVE_FLOAT__COMPLEX",
    "HAVE_DOUBLE__COMPLEX",
    "HAVE_LONG_DOUBLE__COMPLEX",
    "OMPI_MPI_OFFSET_SIZE",
    "OMPI_BUILD_CXX_BINDINGS",
    "OMPI_CXX_SUPPORTS_2D_CONST_CAST",
    "OMPI_PARAM_CHECK",
    "OMPI_HAVE_CXX_EXCEPTION_SUPPORT",
    "OMPI_MAJOR_VERSION",
    "OMPI_MINOR_VERSION",
    "OMPI_RELEASE_VERSION",
    "OPAL_C_HAVE_VISIBILITY",
    "OMPI_PROVIDE_MPI_FILE_INTERFACE",
    "MPI_VERSION",
    "MPI_SUBVERSION",
    //OPENMPI version 1.6.5 :
    "OMPI_WANT_CXX_BINDINGS",
    "OMPI_WANT_F77_BINDINGS",
    "OMPI_WANT_F90_BINDINGS",
    //OPENMPI version 1.4.3
    "OMPI_STDC_HEADERS",
    "OMPI_HAVE_SYSTIME_H",
    "OMPI_HAVE_SYS_SYNCH_H",
    "OMPI_HAVE_LONG_LONG",
    "OMPI_SIZEOF_BOOL",
    "OMPI_SIZEOF_INT"
    ];

immutable types = [
    "OMPI_MPI_OFFSET_TYPE",
    "OMPI_OFFSET_DATATYPE",
    "OMPI_MPI_COUNT_TYPE",
    "OPAL_PTRDIFF_TYPE",
    "OMPI_MPI_AINT_TYPE",
    "ompi_fortran_bogus_type_t",
    "ompi_fortran_integer_t",
    "MPI_Fint",
    //OPENMPI version 1.4.3
    "OMPI_PTRDIFF_TYPE"
];
