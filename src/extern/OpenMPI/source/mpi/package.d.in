/*
 * Copyright (c) 2004-2005 The Trustees of Indiana University and Indiana
 *                         University Research and Technology
 *                         Corporation.  All rights reserved.
 * Copyright (c) 2004-2013 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2004-2007 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 * Copyright (c) 2004-2005 The Regents of the University of California.
 *                         All rights reserved.
 * Copyright (c) 2007-2012 Cisco Systems, Inc.  All rights reserved.
 * Copyright (c) 2008-2009 Sun Microsystems, Inc.  All rights reserved.
 * Copyright (c) 2009-2012 Oak Rigde National Laboratory.  All rights reserved.
 * Copyright (c) 2011      Sandia National Laboratories. All rights reserved.
 * Copyright (c) 2012-2014 Los Alamos Nat Security, LLC. All rights reserved.
 * Copyright (c) 2011-2013 INRIA.  All rights reserved.
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

module mpi;

import core.stdc.config;

extern (C):

// BEGIN AUTO

// END AUTO

static if(OMPI_MAJOR_VERSION == 1 && OMPI_MINOR_VERSION < 8)
{
    pragma(msg, "Version of Open MPI library is too early.");
}
enum OPEN_MPI = 1;

/*
 * To accomodate programs written for MPI implementations that use a
 * straight ROMIO import
 */
alias MPIO_Request = MPI_Request;
alias MPIO_Test = MPI_Test;
alias MPIO_Wait = MPI_Wait;

/*
 * When initializing global pointers to Open MPI internally-defined
 * structs, some compilers warn about type-punning to incomplete
 * types.  Therefore, when full struct definitions are unavailable
 * (when not building Open MPI), cast to an opaque (void *) pointer to
 * disable any strict-aliasing optimizations.  Don't cast to (void *)
 * when building Open MPI so that we actually get the benefit of type
 * checking (because we *do* have the full type definitions available
 * when building OMPI).
 */
//#if !OMPI_BUILDING
//#define OMPI_PREDEFINED_GLOBAL(type, global) ((type) ((void *) &(global)))
//#else
//#define OMPI_PREDEFINED_GLOBAL(type, global) ((type) &(global))
//#endif
To OMPI_PREDEFINED_GLOBAL(To, alias global)()
{
    return cast(To)&global;
}

/*
 * MPI_Status
 */
struct ompi_status_public_t {
    /* These fields are publicly defined in the MPI specification.
       User applications may freely read from these fields. */
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;
    /* The following two fields are internal to the Open MPI
       implementation and should not be accessed by MPI applications.
       They are subject to change at any time.  These are not the
       droids you're looking for. */
    private
    {
        int _cancelled;
        size_t _ucount;
    }
}

/*
 * Typedefs (aliases)
 */
static if(OMPI_MAJOR_VERSION == 1)
    alias MPI_Aint = OPAL_PTRDIFF_TYPE;
else
    alias MPI_Aint = OMPI_MPI_AINT_TYPE;
alias MPI_Offset = OMPI_MPI_OFFSET_TYPE;
alias MPI_Count = OMPI_MPI_COUNT_TYPE;
struct ompi_communicator_t;
struct ompi_datatype_t;
struct ompi_errhandler_t;
struct ompi_file_t;
struct ompi_group_t;
struct ompi_info_t;
struct ompi_op_t;
struct ompi_request_t;
struct ompi_message_t;
// struct ompi_status_public_t;
struct ompi_win_t;
struct mca_base_var_enum_t;
struct ompi_mpit_cvar_handle_t;
struct mca_base_pvar_handle_t;
struct mca_base_pvar_session_t;

alias MPI_Comm = ompi_communicator_t*;
alias MPI_Datatype = ompi_datatype_t*;
alias MPI_Errhandler = ompi_errhandler_t*;
alias MPI_File = ompi_file_t*;
alias MPI_Group = ompi_group_t*;
alias MPI_Info = ompi_info_t*;
alias MPI_Op = ompi_op_t*;
alias MPI_Request = ompi_request_t*;
alias MPI_Message = ompi_message_t*;
alias MPI_Status = ompi_status_public_t;
alias MPI_Win = ompi_win_t*;
alias MPI_T_enum = mca_base_var_enum_t*;
alias MPI_T_cvar_handle = ompi_mpit_cvar_handle_t*;
alias MPI_T_pvar_handle = mca_base_pvar_handle_t*;
alias MPI_T_pvar_session = mca_base_pvar_session_t*;

/*
 * User typedefs
 */
alias MPI_Copy_function = int function(MPI_Comm, int, void *, void *, void *, int *);

alias MPI_Delete_function = int function(MPI_Comm, int, void *, void *);
alias MPI_Datarep_extent_function = int function(MPI_Datatype, MPI_Aint *, void *);
alias MPI_Datarep_conversion_function = int function(void *, MPI_Datatype, int, void *, MPI_Offset, void *);
alias MPI_Comm_errhandler_function = void function(MPI_Comm *, int *, ...);
deprecated alias MPI_Comm_errhandler_fn = MPI_Comm_errhandler_function;


/* This is a little hackish, but errhandler.h needs space for a
   MPI_File_errhandler_fn.  While it could just be removed, this
   allows us to maintain a stable ABI within OMPI, at least for
   apps that don't use MPI I/O. */
alias ompi_file_errhandler_fn = void function(MPI_File *, int *, ...);
alias MPI_File_errhandler_function = ompi_file_errhandler_fn;
alias MPI_Win_errhandler_function = void function(MPI_Win*, int*, ...);
deprecated alias MPI_Win_errhandler_fn = MPI_Win_errhandler_function;
alias MPI_Handler_function = void function(MPI_Comm *, int *, ...);
alias MPI_User_function = void function(void*, void*, int*, MPI_Datatype*);
alias MPI_Comm_copy_attr_function = int function(MPI_Comm, int, void *, void *, void *, int *);

alias MPI_Comm_delete_attr_function = int function(MPI_Comm, int, void *, void *);
alias MPI_Type_copy_attr_function = int function(MPI_Datatype, int, void *, void *, void *, int *);

alias MPI_Type_delete_attr_function = int function(MPI_Datatype, int, void *, void *);

alias MPI_Win_copy_attr_function = int function(MPI_Win, int, void *, void *, void *, int *);

alias MPI_Win_delete_attr_function = int function(MPI_Win, int, void *, void *);
alias MPI_Grequest_query_function = int function(void *, MPI_Status *);
alias MPI_Grequest_free_function = int function(void *);
alias MPI_Grequest_cancel_function = int function(void *, int);

/*
 * Miscellaneous constants
 */
enum MPI_ANY_SOURCE         = -1;                      /* match any source rank */
enum MPI_PROC_NULL          = -2;                      /* rank of null process */
enum MPI_ROOT               = -4;                      /* special value for intercomms */
enum MPI_ANY_TAG            = -1;                      /* match any message tag */
enum MPI_MAX_PROCESSOR_NAME = OPAL_MAX_PROCESSOR_NAME; /* max proc. name length */
enum MPI_MAX_ERROR_STRING   = OPAL_MAX_ERROR_STRING;   /* max error message length */
enum MPI_MAX_OBJECT_NAME    = OPAL_MAX_OBJECT_NAME;    /* max object name length */
enum MPI_MAX_LIBRARY_VERSION_STRING = 256;             /* max length of library version string */
enum MPI_UNDEFINED          = -32766;                  /* undefined stuff */
enum MPI_DIST_GRAPH         = 3;
enum MPI_CART               = 1;                       /* cartesian topology */
enum MPI_GRAPH              = 2;                       /* graph topology */
enum MPI_KEYVAL_INVALID     = -1;                      /* invalid key value */

/*
 * More constants
 */
enum MPI_UNWEIGHTED           = cast(void*) 2;      /* unweighted graph */
enum MPI_WEIGHTS_EMPTY        = cast(void*) 3;      /* empty weights */
enum MPI_BOTTOM               = cast(void*) 0;      /* base reference address */
enum MPI_IN_PLACE             = cast(void*) 1;      /* in place buffer */
enum MPI_BSEND_OVERHEAD       = 128;                /* size of bsend header + ptr */
enum MPI_MAX_INFO_KEY         = OPAL_MAX_INFO_KEY;  /* max info key length */
enum MPI_MAX_INFO_VAL         = OPAL_MAX_INFO_VAL;  /* max info value length */
enum MPI_ARGV_NULL            = cast(char**) 0;     /* NULL argument vector */
enum MPI_ARGVS_NULL           = cast(char***) 0;    /* NULL argument vectors */
enum MPI_ERRCODES_IGNORE      = cast(int*) 0;       /* don't return error codes */
enum MPI_MAX_PORT_NAME        = OPAL_MAX_PORT_NAME; /* max port name length */
enum MPI_ORDER_C              = 0;                  /* C row major order */
enum MPI_ORDER_FORTRAN        = 1;                  /* Fortran column major order */
enum MPI_DISTRIBUTE_BLOCK     = 0;                  /* block distribution */
enum MPI_DISTRIBUTE_CYCLIC    = 1;                  /* cyclic distribution */
enum MPI_DISTRIBUTE_NONE      = 2;                  /* not distributed */
enum MPI_DISTRIBUTE_DFLT_DARG = -1;                 /* default distribution arg */

/*
 * Since these values are arbitrary to Open MPI; we might as well make
 * them the same as ROMIO for ease of mapping.  These values taken
 * from ROMIO's mpio.h file.
 */
enum MPI_MODE_CREATE          =   1; /* ADIO_CREATE */
enum MPI_MODE_RDONLY          =   2; /* ADIO_RDONLY */
enum MPI_MODE_WRONLY          =   4; /* ADIO_WRONLY  */
enum MPI_MODE_RDWR            =   8; /* ADIO_RDWR  */
enum MPI_MODE_DELETE_ON_CLOSE =  16; /* ADIO_DELETE_ON_CLOSE */
enum MPI_MODE_UNIQUE_OPEN     =  32; /* ADIO_UNIQUE_OPEN */
enum MPI_MODE_EXCL            =  64; /* ADIO_EXCL */
enum MPI_MODE_APPEND          = 128; /* ADIO_APPEND */
enum MPI_MODE_SEQUENTIAL      = 256; /* ADIO_SEQUENTIAL */

enum MPI_DISPLACEMENT_CURRENT = -54278278;

enum MPI_SEEK_SET             = 600;
enum MPI_SEEK_CUR             = 602;
enum MPI_SEEK_END             = 604;

/* Max data representation length */
enum MPI_MAX_DATAREP_STRING   = OPAL_MAX_DATAREP_STRING;

/*
 * MPI-2 One-Sided Communications asserts
 */
enum MPI_MODE_NOCHECK        =    1;
enum MPI_MODE_NOPRECEDE      =    2;
enum MPI_MODE_NOPUT          =    4;
enum MPI_MODE_NOSTORE        =    8;
enum MPI_MODE_NOSUCCEED      =   16;

enum MPI_LOCK_EXCLUSIVE      =    1;
enum MPI_LOCK_SHARED         =    2;

enum MPI_WIN_FLAVOR_CREATE   =    1;
enum MPI_WIN_FLAVOR_ALLOCATE =    2;
enum MPI_WIN_FLAVOR_DYNAMIC  =    3;
enum MPI_WIN_FLAVOR_SHARED   =    4;

enum MPI_WIN_UNIFIED         =    0;
enum MPI_WIN_SEPARATE        =    1;

/*
 * Predefined attribute keyvals
 *
 * DO NOT CHANGE THE ORDER WITHOUT ALSO CHANGING THE ORDER IN
 * src/attribute/attribute_predefined.c and mpif.h.in.
 */
enum {
    /* MPI-1 */
    MPI_TAG_UB,
    MPI_HOST,
    MPI_IO,
    MPI_WTIME_IS_GLOBAL,

    /* MPI-2 */
    MPI_APPNUM,
    MPI_LASTUSEDCODE,
    MPI_UNIVERSE_SIZE,
    MPI_WIN_BASE,
    MPI_WIN_SIZE,
    MPI_WIN_DISP_UNIT,

    MPI_WIN_CREATE_FLAVOR,
    MPI_WIN_MODEL,

    /* Even though these four are IMPI attributes, they need to be there
       for all MPI jobs */
    IMPI_CLIENT_SIZE,
    IMPI_CLIENT_COLOR,
    IMPI_HOST_SIZE,
    IMPI_HOST_COLOR
}

/*
 * Error classes and codes
 * Do not change the values of these without also modifying mpif.h.in.
 */
enum MPI_SUCCESS                   =  0;
enum MPI_ERR_BUFFER                =  1;
enum MPI_ERR_COUNT                 =  2;
enum MPI_ERR_TYPE                  =  3;
enum MPI_ERR_TAG                   =  4;
enum MPI_ERR_COMM                  =  5;
enum MPI_ERR_RANK                  =  6;
enum MPI_ERR_REQUEST               =  7;
enum MPI_ERR_ROOT                  =  8;
enum MPI_ERR_GROUP                 =  9;
enum MPI_ERR_OP                    = 10;
enum MPI_ERR_TOPOLOGY              = 11;
enum MPI_ERR_DIMS                  = 12;
enum MPI_ERR_ARG                   = 13;
enum MPI_ERR_UNKNOWN               = 14;
enum MPI_ERR_TRUNCATE              = 15;
enum MPI_ERR_OTHER                 = 16;
enum MPI_ERR_INTERN                = 17;
enum MPI_ERR_IN_STATUS             = 18;
enum MPI_ERR_PENDING               = 19;
enum MPI_ERR_ACCESS                = 20;
enum MPI_ERR_AMODE                 = 21;
enum MPI_ERR_ASSERT                = 22;
enum MPI_ERR_BAD_FILE              = 23;
enum MPI_ERR_BASE                  = 24;
enum MPI_ERR_CONVERSION            = 25;
enum MPI_ERR_DISP                  = 26;
enum MPI_ERR_DUP_DATAREP           = 27;
enum MPI_ERR_FILE_EXISTS           = 28;
enum MPI_ERR_FILE_IN_USE           = 29;
enum MPI_ERR_FILE                  = 30;
enum MPI_ERR_INFO_KEY              = 31;
enum MPI_ERR_INFO_NOKEY            = 32;
enum MPI_ERR_INFO_VALUE            = 33;
enum MPI_ERR_INFO                  = 34;
enum MPI_ERR_IO                    = 35;
enum MPI_ERR_KEYVAL                = 36;
enum MPI_ERR_LOCKTYPE              = 37;
enum MPI_ERR_NAME                  = 38;
enum MPI_ERR_NO_MEM                = 39;
enum MPI_ERR_NOT_SAME              = 40;
enum MPI_ERR_NO_SPACE              = 41;
enum MPI_ERR_NO_SUCH_FILE          = 42;
enum MPI_ERR_PORT                  = 43;
enum MPI_ERR_QUOTA                 = 44;
enum MPI_ERR_READ_ONLY             = 45;
enum MPI_ERR_RMA_CONFLICT          = 46;
enum MPI_ERR_RMA_SYNC              = 47;
enum MPI_ERR_SERVICE               = 48;
enum MPI_ERR_SIZE                  = 49;
enum MPI_ERR_SPAWN                 = 50;
enum MPI_ERR_UNSUPPORTED_DATAREP   = 51;
enum MPI_ERR_UNSUPPORTED_OPERATION = 52;
enum MPI_ERR_WIN                   = 53;
enum MPI_T_ERR_MEMORY              = 54;
enum MPI_T_ERR_NOT_INITIALIZED     = 55;
enum MPI_T_ERR_CANNOT_INIT         = 56;
enum MPI_T_ERR_INVALID_INDEX       = 57;
enum MPI_T_ERR_INVALID_ITEM        = 58;
enum MPI_T_ERR_INVALID_HANDLE      = 59;
enum MPI_T_ERR_OUT_OF_HANDLES      = 60;
enum MPI_T_ERR_OUT_OF_SESSIONS     = 61;
enum MPI_T_ERR_INVALID_SESSION     = 62;
enum MPI_T_ERR_CVAR_SET_NOT_NOW    = 63;
enum MPI_T_ERR_CVAR_SET_NEVER      = 64;
enum MPI_T_ERR_PVAR_NO_STARTSTOP   = 65;
enum MPI_T_ERR_PVAR_NO_WRITE       = 66;
enum MPI_T_ERR_PVAR_NO_ATOMIC      = 67;
enum MPI_ERR_RMA_RANGE             = 68;
enum MPI_ERR_RMA_ATTACH            = 69;
enum MPI_ERR_RMA_FLAVOR            = 70;
enum MPI_ERR_RMA_SHARED            = 71;
enum MPI_T_ERR_INVALID             = 72;
enum MPI_T_ERR_INVALID_NAME        = 73;
    
/* Per MPI-3 p349 47, MPI_ERR_LASTCODE must be >= the last predefined
   MPI_ERR_<foo> code.  Set the last code to allow some room for adding
   error codes without breaking ABI. */
enum MPI_ERR_LASTCODE              = 92;

/*
 * Comparison results.  Don't change the order of these, the group
 * comparison functions rely on it.
 * Do not change the order of these without also modifying mpif.h.in.
 */
enum
{
    MPI_IDENT,
    MPI_CONGRUENT,
    MPI_SIMILAR,
    MPI_UNEQUAL
}

/*
 * MPI_Init_thread constants
 * Do not change the order of these without also modifying mpif.h.in.
 */
enum
{
    MPI_THREAD_SINGLE,
    MPI_THREAD_FUNNELED,
    MPI_THREAD_SERIALIZED,
    MPI_THREAD_MULTIPLE
}

/*
 * Datatype combiners.
 * Do not change the order of these without also modifying mpif.h.in.
 */
enum
{
    MPI_COMBINER_NAMED,
    MPI_COMBINER_DUP,
    MPI_COMBINER_CONTIGUOUS,
    MPI_COMBINER_VECTOR,
    MPI_COMBINER_HVECTOR_INTEGER,
    MPI_COMBINER_HVECTOR,
    MPI_COMBINER_INDEXED,
    MPI_COMBINER_HINDEXED_INTEGER,
    MPI_COMBINER_HINDEXED,
    MPI_COMBINER_INDEXED_BLOCK,
    MPI_COMBINER_STRUCT_INTEGER,
    MPI_COMBINER_STRUCT,
    MPI_COMBINER_SUBARRAY,
    MPI_COMBINER_DARRAY,
    MPI_COMBINER_F90_REAL,
    MPI_COMBINER_F90_COMPLEX,
    MPI_COMBINER_F90_INTEGER,
    MPI_COMBINER_RESIZED,
    MPI_COMBINER_HINDEXED_BLOCK
}

/*
 * Communicator split type constants.
 * Do not change the order of these without also modifying mpif.h.in
 * (see also mpif-common.h.fin).
 */
enum {
      MPI_COMM_TYPE_SHARED,
      OMPI_COMM_TYPE_HWTHREAD,
      OMPI_COMM_TYPE_CORE,
      OMPI_COMM_TYPE_L1CACHE,
      OMPI_COMM_TYPE_L2CACHE,
      OMPI_COMM_TYPE_L3CACHE,
      OMPI_COMM_TYPE_SOCKET,
      OMPI_COMM_TYPE_NUMA,
      OMPI_COMM_TYPE_BOARD,
      OMPI_COMM_TYPE_HOST,
      OMPI_COMM_TYPE_CU,
      OMPI_COMM_TYPE_CLUSTER
}
alias OMPI_COMM_TYPE_NODE = MPI_COMM_TYPE_SHARED;
    
/*
 * MPIT Verbosity Levels
 */
enum {
      MPI_T_VERBOSITY_USER_BASIC,
      MPI_T_VERBOSITY_USER_DETAIL,
      MPI_T_VERBOSITY_USER_ALL,
      MPI_T_VERBOSITY_TUNER_BASIC,
      MPI_T_VERBOSITY_TUNER_DETAIL,
      MPI_T_VERBOSITY_TUNER_ALL,
      MPI_T_VERBOSITY_MPIDEV_BASIC,
      MPI_T_VERBOSITY_MPIDEV_DETAIL,
      MPI_T_VERBOSITY_MPIDEV_ALL
}
    
/*
 * MPIT Scopes
 */
enum {
      MPI_T_SCOPE_CONSTANT,
      MPI_T_SCOPE_READONLY,
      MPI_T_SCOPE_LOCAL,
      MPI_T_SCOPE_GROUP,
      MPI_T_SCOPE_GROUP_EQ,
      MPI_T_SCOPE_ALL,
      MPI_T_SCOPE_ALL_EQ
}
    
/*
 * MPIT Object Binding
 */
enum {
      MPI_T_BIND_NO_OBJECT,
      MPI_T_BIND_MPI_COMM,
      MPI_T_BIND_MPI_DATATYPE,
      MPI_T_BIND_MPI_ERRHANDLER,
      MPI_T_BIND_MPI_FILE,
      MPI_T_BIND_MPI_GROUP,
      MPI_T_BIND_MPI_OP,
      MPI_T_BIND_MPI_REQUEST,
      MPI_T_BIND_MPI_WIN,
      MPI_T_BIND_MPI_MESSAGE,
      MPI_T_BIND_MPI_INFO
}
    
/*
 * MPIT pvar classes
 */
enum {
      MPI_T_PVAR_CLASS_STATE,
      MPI_T_PVAR_CLASS_LEVEL,
      MPI_T_PVAR_CLASS_SIZE,
      MPI_T_PVAR_CLASS_PERCENTAGE,
      MPI_T_PVAR_CLASS_HIGHWATERMARK,
      MPI_T_PVAR_CLASS_LOWWATERMARK,
      MPI_T_PVAR_CLASS_COUNTER,
      MPI_T_PVAR_CLASS_AGGREGATE,
      MPI_T_PVAR_CLASS_TIMER,
      MPI_T_PVAR_CLASS_GENERIC
}

/*
 * NULL handles
 */
alias MPI_GROUP_NULL = OMPI_PREDEFINED_GLOBAL!(MPI_Group, ompi_mpi_group_null);
alias MPI_COMM_NULL = OMPI_PREDEFINED_GLOBAL!(MPI_Comm, ompi_mpi_comm_null);
alias MPI_REQUEST_NULL = OMPI_PREDEFINED_GLOBAL!(MPI_Request, ompi_request_null);
static if(OMPI_MAJOR_VERSION == 1 && OMPI_MINOR_VERSION >= 8)
    alias MPI_MESSAGE_NULL = OMPI_PREDEFINED_GLOBAL!(MPI_Message, ompi_message_null);
alias MPI_OP_NULL = OMPI_PREDEFINED_GLOBAL!(MPI_Op, ompi_mpi_op_null);
alias MPI_ERRHANDLER_NULL = OMPI_PREDEFINED_GLOBAL!(MPI_Errhandler, ompi_mpi_errhandler_null);
alias MPI_INFO_NULL = OMPI_PREDEFINED_GLOBAL!(MPI_Info, ompi_mpi_info_null);
alias MPI_WIN_NULL = OMPI_PREDEFINED_GLOBAL!(MPI_Win, ompi_mpi_win_null);
alias MPI_FILE_NULL = OMPI_PREDEFINED_GLOBAL!(MPI_File, ompi_mpi_file_null);
enum MPI_T_ENUM_NULL = cast(MPI_T_enum) null;
    
/*
 * MPI_INFO_ENV handle
 */
alias MPI_INFO_ENV = OMPI_PREDEFINED_GLOBAL!(MPI_Info, ompi_mpi_info_env);

enum MPI_STATUS_IGNORE = cast(MPI_Status*) 0;
enum MPI_STATUSES_IGNORE = cast(MPI_Status*) 0;

/*
 * Special MPI_T handles
 */
enum MPI_T_PVAR_ALL_HANDLES = cast(MPI_T_pvar_handle) -1;
enum MPI_T_PVAR_HANDLE_NULL = cast(MPI_T_pvar_handle) 0;
enum MPI_T_CVAR_HANDLE_NULL = cast(MPI_T_cvar_handle) 0;

/* MPI-2 specifies that the name "MPI_TYPE_NULL_DELETE_FN" (and all
   related friends) must be accessible in C, C++, and Fortran. This is
   unworkable if the back-end Fortran compiler uses all caps for its
   linker symbol convention -- it results in two functions with
   different signatures that have the same name (i.e., both C and
   Fortran use the symbol MPI_TYPE_NULL_DELETE_FN).  So we have to
   #define the C names to be something else, so that they names are
   *accessed* through MPI_TYPE_NULL_DELETE_FN, but their actual symbol
   name is different.

   However, this file is included when the fortran wrapper functions
   are compiled in Open MPI, so we do *not* want these #defines in
   this case (i.e., we need the Fortran wrapper function to be
   compiled as MPI_TYPE_NULL_DELETE_FN).  So add some #if kinds of
   protection for this case. */

static if (!is(typeof(OMPI_COMPILING_FORTRAN_WRAPPERS)))
{
    alias MPI_NULL_DELETE_FN = OMPI_C_MPI_NULL_DELETE_FN;
    alias MPI_NULL_COPY_FN = OMPI_C_MPI_NULL_COPY_FN;
    alias MPI_DUP_FN = OMPI_C_MPI_DUP_FN;
    
    alias MPI_TYPE_NULL_DELETE_FN = OMPI_C_MPI_TYPE_NULL_DELETE_FN;
    alias MPI_TYPE_NULL_COPY_FN = OMPI_C_MPI_TYPE_NULL_COPY_FN;
    alias MPI_TYPE_DUP_FN = OMPI_C_MPI_TYPE_DUP_FN;
    
    alias MPI_COMM_NULL_DELETE_FN = OMPI_C_MPI_COMM_NULL_DELETE_FN;
    alias MPI_COMM_NULL_COPY_FN = OMPI_C_MPI_COMM_NULL_COPY_FN;
    alias MPI_COMM_DUP_FN = OMPI_C_MPI_COMM_DUP_FN;
    
    alias MPI_WIN_NULL_DELETE_FN = OMPI_C_MPI_WIN_NULL_DELETE_FN;
    alias MPI_WIN_NULL_COPY_FN = OMPI_C_MPI_WIN_NULL_COPY_FN;
    alias MPI_WIN_DUP_FN = OMPI_C_MPI_WIN_DUP_FN;
    
    /* MPI_CONVERSION_FN_NULL is a sentinel value, but it has to be large
       enough to be the same size as a valid function pointer.  It
       therefore shares many characteristics between Fortran constants and
       Fortran sentinel functions.  For example, it shares the problem of
       having Fortran compilers have all-caps versions of the symbols that
       must be able to be present, and therefore has to be in this
       conditional block in mpi.h. */
    enum MPI_CONVERSION_FN_NULL = cast(MPI_Datarep_conversion_function*) 0;
}
else
{
    enum _TEST_ = OMPI_COMPILING_FORTRAN_WRAPPERS;
}

int OMPI_C_MPI_TYPE_NULL_DELETE_FN( MPI_Datatype datatype,
        int type_keyval,
        void* attribute_val_out,
        void* extra_state );
int OMPI_C_MPI_TYPE_NULL_COPY_FN( MPI_Datatype datatype,
        int type_keyval,
        void* extra_state,
        void* attribute_val_in,
        void* attribute_val_out,
        int* flag );
int OMPI_C_MPI_TYPE_DUP_FN( MPI_Datatype datatype,
        int type_keyval,
        void* extra_state,
        void* attribute_val_in,
        void* attribute_val_out,
        int* flag );
int OMPI_C_MPI_COMM_NULL_DELETE_FN( MPI_Comm comm,
        int comm_keyval,
        void* attribute_val_out,
        void* extra_state );
int OMPI_C_MPI_COMM_NULL_COPY_FN( MPI_Comm comm,
        int comm_keyval,
        void* extra_state,
        void* attribute_val_in,
        void* attribute_val_out,
        int* flag );
int OMPI_C_MPI_COMM_DUP_FN( MPI_Comm comm, int comm_keyval,
        void* extra_state,
        void* attribute_val_in,
        void* attribute_val_out,
        int* flag );
int OMPI_C_MPI_NULL_DELETE_FN( MPI_Comm comm, int comm_keyval,
        void* attribute_val_out,
        void* extra_state );
int OMPI_C_MPI_NULL_COPY_FN( MPI_Comm comm, int comm_keyval,
        void* extra_state,
        void* attribute_val_in,
        void* attribute_val_out,
        int* flag );
int OMPI_C_MPI_DUP_FN( MPI_Comm comm, int comm_keyval,
        void* extra_state,
        void* attribute_val_in,
        void* attribute_val_out,
        int* flag );
int OMPI_C_MPI_WIN_NULL_DELETE_FN( MPI_Win window,
        int win_keyval,
        void* attribute_val_out,
        void* extra_state );
int OMPI_C_MPI_WIN_NULL_COPY_FN( MPI_Win window, int win_keyval,
        void* extra_state,
        void* attribute_val_in,
        void* attribute_val_out,
        int* flag );
int OMPI_C_MPI_WIN_DUP_FN( MPI_Win window, int win_keyval,
        void* extra_state,
        void* attribute_val_in,
        void* attribute_val_out,
        int* flag );

/*
 * External variables
 *
 * The below externs use the ompi_predefined_xxx_t structures to maintain
 * back compatibility between MPI library versions.
 * See ompi/communicator/communicator.h comments with struct ompi_communicator_t
 * for full explanation why we chose to use the ompi_predefined_xxx_t structure.
 */
extern __gshared
{
    struct ompi_predefined_communicator_t {}
    ompi_predefined_communicator_t ompi_mpi_comm_world;
    ompi_predefined_communicator_t ompi_mpi_comm_self;
    ompi_predefined_communicator_t ompi_mpi_comm_null;

    struct ompi_predefined_group_t {}
    ompi_predefined_group_t ompi_mpi_group_empty;
    ompi_predefined_group_t ompi_mpi_group_null;

    struct ompi_predefined_request_t {}
    ompi_predefined_request_t ompi_request_null;

    struct ompi_predefined_message_t {}
    ompi_predefined_message_t ompi_message_null;
    ompi_predefined_message_t ompi_message_no_proc;

    struct ompi_predefined_op_t {}
    ompi_predefined_op_t ompi_mpi_op_null;
    ompi_predefined_op_t ompi_mpi_op_min;
    ompi_predefined_op_t ompi_mpi_op_max;
    ompi_predefined_op_t ompi_mpi_op_sum;
    ompi_predefined_op_t ompi_mpi_op_prod;
    ompi_predefined_op_t ompi_mpi_op_land;
    ompi_predefined_op_t ompi_mpi_op_band;
    ompi_predefined_op_t ompi_mpi_op_lor;
    ompi_predefined_op_t ompi_mpi_op_bor;
    ompi_predefined_op_t ompi_mpi_op_lxor;
    ompi_predefined_op_t ompi_mpi_op_bxor;
    ompi_predefined_op_t ompi_mpi_op_maxloc;
    ompi_predefined_op_t ompi_mpi_op_minloc;
    ompi_predefined_op_t ompi_mpi_op_replace;
    ompi_predefined_op_t ompi_mpi_op_no_op;

    struct ompi_predefined_datatype_t {}
    ompi_predefined_datatype_t ompi_mpi_datatype_null;

    // ompi_predefined_datatype_t ompi_mpi_lb; // not in recent MPI
    // ompi_predefined_datatype_t ompi_mpi_ub; // not in recent MPI
    ompi_predefined_datatype_t ompi_mpi_char;
    ompi_predefined_datatype_t ompi_mpi_signed_char;
    ompi_predefined_datatype_t ompi_mpi_unsigned_char;
    ompi_predefined_datatype_t ompi_mpi_byte;
    ompi_predefined_datatype_t ompi_mpi_short;
    ompi_predefined_datatype_t ompi_mpi_unsigned_short;
    ompi_predefined_datatype_t ompi_mpi_int;
    ompi_predefined_datatype_t ompi_mpi_unsigned;
    ompi_predefined_datatype_t ompi_mpi_long;
    ompi_predefined_datatype_t ompi_mpi_unsigned_long;
    ompi_predefined_datatype_t ompi_mpi_long_long_int;
    ompi_predefined_datatype_t ompi_mpi_unsigned_long_long;
    ompi_predefined_datatype_t ompi_mpi_float;
    ompi_predefined_datatype_t ompi_mpi_double;
    ompi_predefined_datatype_t ompi_mpi_long_double;
    ompi_predefined_datatype_t ompi_mpi_wchar;
    ompi_predefined_datatype_t ompi_mpi_packed;

    /*
     * Following are the C++/C99 datatypes
     */
    ompi_predefined_datatype_t ompi_mpi_cxx_bool;
    ompi_predefined_datatype_t ompi_mpi_cxx_cplex;
    ompi_predefined_datatype_t ompi_mpi_cxx_dblcplex;
    ompi_predefined_datatype_t ompi_mpi_cxx_ldblcplex;

    /*
     * Following are the Fortran datatypes
     */
    ompi_predefined_datatype_t ompi_mpi_logical;
    ompi_predefined_datatype_t ompi_mpi_character;
    ompi_predefined_datatype_t ompi_mpi_integer;
    ompi_predefined_datatype_t ompi_mpi_real;
    ompi_predefined_datatype_t ompi_mpi_dblprec;
    ompi_predefined_datatype_t ompi_mpi_cplex;
    ompi_predefined_datatype_t ompi_mpi_dblcplex;
    ompi_predefined_datatype_t ompi_mpi_ldblcplex;

    /* Aggregate struct datatypes are not const */
    ompi_predefined_datatype_t ompi_mpi_2int;
    ompi_predefined_datatype_t ompi_mpi_2integer;
    ompi_predefined_datatype_t ompi_mpi_2real;
    ompi_predefined_datatype_t ompi_mpi_2dblprec;
    ompi_predefined_datatype_t ompi_mpi_2cplex;
    ompi_predefined_datatype_t ompi_mpi_2dblcplex;

    ompi_predefined_datatype_t ompi_mpi_float_int;
    ompi_predefined_datatype_t ompi_mpi_double_int;
    ompi_predefined_datatype_t ompi_mpi_longdbl_int;
    ompi_predefined_datatype_t ompi_mpi_short_int;
    ompi_predefined_datatype_t ompi_mpi_long_int;

    //OpenMPI used to version these based on OMPI_HAVE_FORTRAN_LOGICALx,
    //but it doesn't matter, the actual MPI_LOGICALx symbols are versioned
    //and nothing is emitted for unused extern variables
    /* Optional MPI2 datatypes, always declared and defined, but not "exported" as MPI_LOGICAL1 */
    ompi_predefined_datatype_t ompi_mpi_logical1;
    ompi_predefined_datatype_t ompi_mpi_logical2;
    ompi_predefined_datatype_t ompi_mpi_logical4;
    ompi_predefined_datatype_t ompi_mpi_logical8;
    ompi_predefined_datatype_t ompi_mpi_integer1;
    ompi_predefined_datatype_t ompi_mpi_integer2;
    ompi_predefined_datatype_t ompi_mpi_integer4;
    ompi_predefined_datatype_t ompi_mpi_integer8;
    ompi_predefined_datatype_t ompi_mpi_integer16;
    ompi_predefined_datatype_t ompi_mpi_real2;
    ompi_predefined_datatype_t ompi_mpi_real4;
    ompi_predefined_datatype_t ompi_mpi_real8;
    ompi_predefined_datatype_t ompi_mpi_real16;
    ompi_predefined_datatype_t ompi_mpi_complex8;
    ompi_predefined_datatype_t ompi_mpi_complex16;
    ompi_predefined_datatype_t ompi_mpi_complex32;

    /* New datatypes from the MPI 2.2 standard */
    ompi_predefined_datatype_t ompi_mpi_int8_t;
    ompi_predefined_datatype_t ompi_mpi_uint8_t;
    ompi_predefined_datatype_t ompi_mpi_int16_t;
    ompi_predefined_datatype_t ompi_mpi_uint16_t;
    ompi_predefined_datatype_t ompi_mpi_int32_t;
    ompi_predefined_datatype_t ompi_mpi_uint32_t;
    ompi_predefined_datatype_t ompi_mpi_int64_t;
    ompi_predefined_datatype_t ompi_mpi_uint64_t;
    ompi_predefined_datatype_t ompi_mpi_aint;
    ompi_predefined_datatype_t ompi_mpi_offset;
    ompi_predefined_datatype_t ompi_mpi_count;
    ompi_predefined_datatype_t ompi_mpi_c_bool;
    ompi_predefined_datatype_t ompi_mpi_c_complex; // not in v4.0.x header
    ompi_predefined_datatype_t ompi_mpi_c_float_complex;
    ompi_predefined_datatype_t ompi_mpi_c_double_complex;
    ompi_predefined_datatype_t ompi_mpi_c_long_double_complex;

    struct ompi_predefined_errhandler_t {}
    ompi_predefined_errhandler_t ompi_mpi_errhandler_null;
    ompi_predefined_errhandler_t ompi_mpi_errors_are_fatal;
    ompi_predefined_errhandler_t ompi_mpi_errors_return;

    struct ompi_predefined_win_t {}
    ompi_predefined_win_t ompi_mpi_win_null;

    struct ompi_predefined_file_t {}
    ompi_predefined_file_t ompi_mpi_file_null;

    struct ompi_predefined_info_t {}
    ompi_predefined_info_t ompi_mpi_info_null;
    ompi_predefined_info_t ompi_mpi_info_env;

    MPI_Fint* MPI_F_STATUS_IGNORE;
    MPI_Fint* MPI_F_STATUSES_IGNORE;
}

/*
 * MPI predefined handles
 */
alias MPI_COMM_WORLD = OMPI_PREDEFINED_GLOBAL!(MPI_Comm, ompi_mpi_comm_world);
alias MPI_COMM_SELF = OMPI_PREDEFINED_GLOBAL!(MPI_Comm, ompi_mpi_comm_self);

alias MPI_GROUP_EMPTY = OMPI_PREDEFINED_GLOBAL!(MPI_Group, ompi_mpi_group_empty);

alias MPI_MESSAGE_NO_PROC = OMPI_PREDEFINED_GLOBAL!(MPI_Message, ompi_message_no_proc);

alias MPI_MAX = OMPI_PREDEFINED_GLOBAL!(MPI_Op, ompi_mpi_op_max);
alias MPI_MIN = OMPI_PREDEFINED_GLOBAL!(MPI_Op, ompi_mpi_op_min);
alias MPI_SUM = OMPI_PREDEFINED_GLOBAL!(MPI_Op, ompi_mpi_op_sum);
alias MPI_PROD = OMPI_PREDEFINED_GLOBAL!(MPI_Op, ompi_mpi_op_prod);
alias MPI_LAND = OMPI_PREDEFINED_GLOBAL!(MPI_Op, ompi_mpi_op_land);
alias MPI_BAND = OMPI_PREDEFINED_GLOBAL!(MPI_Op, ompi_mpi_op_band);
alias MPI_LOR = OMPI_PREDEFINED_GLOBAL!(MPI_Op, ompi_mpi_op_lor);
alias MPI_BOR = OMPI_PREDEFINED_GLOBAL!(MPI_Op, ompi_mpi_op_bor);
alias MPI_LXOR = OMPI_PREDEFINED_GLOBAL!(MPI_Op, ompi_mpi_op_lxor);
alias MPI_BXOR = OMPI_PREDEFINED_GLOBAL!(MPI_Op, ompi_mpi_op_bxor);
alias MPI_MAXLOC = OMPI_PREDEFINED_GLOBAL!(MPI_Op, ompi_mpi_op_maxloc);
alias MPI_MINLOC = OMPI_PREDEFINED_GLOBAL!(MPI_Op, ompi_mpi_op_minloc);
alias MPI_REPLACE = OMPI_PREDEFINED_GLOBAL!(MPI_Op, ompi_mpi_op_replace);
alias MPI_NO_OP = OMPI_PREDEFINED_GLOBAL!(MPI_Op, ompi_mpi_op_no_op);

/* C datatypes */
alias MPI_DATATYPE_NULL = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_datatype_null);
alias MPI_BYTE = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_byte);
alias MPI_PACKED = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_packed);
alias MPI_CHAR = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_char);
alias MPI_SHORT = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_short);
alias MPI_INT = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_int);
alias MPI_LONG = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_long);
alias MPI_FLOAT = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_float);
alias MPI_DOUBLE = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_double);
alias MPI_LONG_DOUBLE = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_long_double);
alias MPI_UNSIGNED_CHAR = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_unsigned_char);
alias MPI_SIGNED_CHAR = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_signed_char);
alias MPI_UNSIGNED_SHORT = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_unsigned_short);
alias MPI_UNSIGNED_LONG = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_unsigned_long);
alias MPI_UNSIGNED = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_unsigned);
alias MPI_FLOAT_INT = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_float_int);
alias MPI_DOUBLE_INT = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_double_int);
alias MPI_LONG_DOUBLE_INT = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_longdbl_int);
alias MPI_LONG_INT = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_long_int);
alias MPI_SHORT_INT = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_short_int);
alias MPI_2INT = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_2int);
// alias MPI_UB = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_ub); // not in recent MPI
// alias MPI_LB = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_lb); // not in recent MPI
alias MPI_WCHAR = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_wchar);
static if (OPAL_HAVE_LONG_LONG)
    private enum _HAVE_LONG = 1;
else
    private enum _HAVE_LONG = 0;
static if (_HAVE_LONG)
{
    alias MPI_LONG_LONG_INT = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_long_long_int);
    alias MPI_LONG_LONG = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_long_long_int);
    alias MPI_UNSIGNED_LONG_LONG = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_unsigned_long_long);
}
// BEGIN AUTO

// END AUTO
alias MPI_2COMPLEX = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_2cplex);
alias MPI_2DOUBLE_COMPLEX = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_2dblcplex);

/* Fortran datatype bindings */
alias MPI_CHARACTER = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_character);
static if(OMPI_MAJOR_VERSION == 1 && OMPI_MINOR_VERSION < 6)
    alias MPI_LOGICAL = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_logic);
else
    alias MPI_LOGICAL = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_logical);
static if (OMPI_HAVE_FORTRAN_LOGICAL1)
    alias MPI_LOGICAL1 = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_logical1);
static if (OMPI_HAVE_FORTRAN_LOGICAL2)
    alias MPI_LOGICAL2 = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_logical2);
static if (OMPI_HAVE_FORTRAN_LOGICAL4)
    alias MPI_LOGICAL4 = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_logical4);
static if (OMPI_HAVE_FORTRAN_LOGICAL8)
    alias MPI_LOGICAL8 = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_logical8);
alias MPI_INTEGER = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_integer);
static if (OMPI_HAVE_FORTRAN_INTEGER1)
    alias MPI_INTEGER1 = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_integer1);
static if (OMPI_HAVE_FORTRAN_INTEGER2)
    alias MPI_INTEGER2 = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_integer2);
static if (OMPI_HAVE_FORTRAN_INTEGER4)
    alias MPI_INTEGER4 = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_integer4);
static if (OMPI_HAVE_FORTRAN_INTEGER8)
    alias MPI_INTEGER8 = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_integer8);
static if (OMPI_HAVE_FORTRAN_INTEGER16)
    alias MPI_INTEGER16 = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_integer16);
alias MPI_REAL = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_real);
static if (OMPI_HAVE_FORTRAN_REAL4)
    alias MPI_REAL4 = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_real4);
static if (OMPI_HAVE_FORTRAN_REAL8)
    alias MPI_REAL8 = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_real8);
static if (OMPI_HAVE_FORTRAN_REAL16)
    alias MPI_REAL16 = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_real16);
alias MPI_DOUBLE_PRECISION = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_dblprec);
alias MPI_COMPLEX = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_cplex);
static if (OMPI_HAVE_FORTRAN_REAL4)
    alias MPI_COMPLEX8 = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_complex8);
static if (OMPI_HAVE_FORTRAN_REAL8)
    alias MPI_COMPLEX16 = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_complex16);
static if (OMPI_HAVE_FORTRAN_REAL16)
    alias MPI_COMPLEX32 = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_complex32);
alias MPI_DOUBLE_COMPLEX = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_dblcplex);
alias MPI_2REAL = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_2real);
alias MPI_2DOUBLE_PRECISION = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_2dblprec);
alias MPI_2INTEGER = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_2integer);

/* New datatypes from the MPI 2.2 standard */
alias MPI_INT8_T                = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_int8_t);
alias MPI_UINT8_T               = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_uint8_t);
alias MPI_INT16_T               = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_int16_t);
alias MPI_UINT16_T              = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_uint16_t);
alias MPI_INT32_T               = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_int32_t);
alias MPI_UINT32_T              = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_uint32_t);
alias MPI_INT64_T               = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_int64_t);
alias MPI_UINT64_T              = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_uint64_t);
alias MPI_AINT                  = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_aint);
alias MPI_OFFSET                = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_offset);
alias MPI_C_BOOL                = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_c_bool);
static if (HAVE_FLOAT__COMPLEX)
{
    alias MPI_C_COMPLEX             = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_c_complex);
    alias MPI_C_FLOAT_COMPLEX       = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_c_float_complex);
 }
static if (HAVE_DOUBLE__COMPLEX)
    alias MPI_C_DOUBLE_COMPLEX      = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_c_double_complex);
static if (HAVE_LONG_DOUBLE__COMPLEX)
    alias MPI_C_LONG_DOUBLE_COMPLEX = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_c_long_double_complex);
alias MPI_CXX_BOOL              = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_cxx_bool);
alias MPI_CXX_FLOAT_COMPLEX     = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_cxx_cplex);
alias MPI_CXX_DOUBLE_COMPLEX    = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_cxx_dblcplex);
alias MPI_CXX_LONG_DOUBLE_COMPLEX = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_cxx_ldblcplex);
    
/* New datatypes from the 3.0 standard */
alias MPI_COUNT                 = OMPI_PREDEFINED_GLOBAL!(MPI_Datatype, ompi_mpi_count);

alias MPI_ERRORS_ARE_FATAL = OMPI_PREDEFINED_GLOBAL!(MPI_Errhandler, ompi_mpi_errors_are_fatal);
alias MPI_ERRORS_RETURN = OMPI_PREDEFINED_GLOBAL!(MPI_Errhandler, ompi_mpi_errors_return);

/* Typeclass definition for MPI_Type_match_size */
enum MPI_TYPECLASS_INTEGER  = 1;
enum MPI_TYPECLASS_REAL     = 2;
enum MPI_TYPECLASS_COMPLEX  = 3;

/* Aint helper macros (MPI-3.1) */
// #define MPI_Aint_add(base, disp) ((MPI_Aint) ((char *) (base) + (disp)))
// #define MPI_Aint_diff(addr1, addr2) ((MPI_Aint) ((char *) (addr1) - (char *) (addr2)))
// #define PMPI_Aint_add(base, disp) MPI_Aint_add(base, disp)
// #define PMPI_Aint_diff(addr1, addr2) MPI_Aint_diff(addr1, addr2)

/*
 * MPI API
 */
//TODO properly version the API from here down
int MPI_Abort(MPI_Comm comm, int errorcode);
int MPI_Accumulate(void* origin_addr, int origin_count, MPI_Datatype origin_datatype,
        int target_rank, MPI_Aint target_disp, int target_count,
        MPI_Datatype target_datatype, MPI_Op op, MPI_Win win);
int MPI_Add_error_class(int* errorclass);
int MPI_Add_error_code(int errorclass, int* errorcode);
int MPI_Add_error_string(int errorcode, char* string);
int MPI_Allgather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount,
        MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Iallgather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount,
        MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request);
int MPI_Allgatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int* recvcounts,
        int* displs, MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Iallgatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int* recvcounts,
        int* displs, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request);
int MPI_Alloc_mem(MPI_Aint size, MPI_Info info,
        void* baseptr);
int MPI_Allreduce(void* sendbuf, void* recvbuf, int count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPI_Iallreduce(void* sendbuf, void* recvbuf, int count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request* request);
int MPI_Alltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount,
        MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Ialltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount,
        MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request);
int MPI_Alltoallv(void* sendbuf, int* sendcounts, int* sdispls,
        MPI_Datatype sendtype, void* recvbuf, int* recvcounts,
        int* rdispls, MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Ialltoallv(void* sendbuf, int* sendcounts, int* sdispls,
        MPI_Datatype sendtype, void* recvbuf, int* recvcounts,
        int* rdispls, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request);
int MPI_Alltoallw(void* sendbuf, int* sendcounts, int* sdispls, MPI_Datatype* sendtypes,
        void* recvbuf, int* recvcounts, int* rdispls, MPI_Datatype* recvtypes,
        MPI_Comm comm);
int MPI_Ialltoallw(void* sendbuf, int* sendcounts, int* sdispls, MPI_Datatype* sendtypes,
        void* recvbuf, int* recvcounts, int* rdispls, MPI_Datatype* recvtypes,
        MPI_Comm comm, MPI_Request* request);
int MPI_Barrier(MPI_Comm comm);
int MPI_Ibarrier(MPI_Comm comm, MPI_Request* request);
int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype,
        int root, MPI_Comm comm);
int MPI_Bsend(void* buf, int count, MPI_Datatype datatype,
        int dest, int tag, MPI_Comm comm);
int MPI_Ibcast(void* buffer, int count, MPI_Datatype datatype,
        int root, MPI_Comm comm,
        MPI_Request* request);
int MPI_Bsend_init(void* buf, int count, MPI_Datatype datatype,
        int dest, int tag, MPI_Comm comm, MPI_Request* request);
int MPI_Buffer_attach(void* buffer, int size);
int MPI_Buffer_detach(void* buffer, int* size);
int MPI_Cancel(MPI_Request* request);
int MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int* coords);
int MPI_Cart_create(MPI_Comm old_comm, int ndims, int* dims,
        int* periods, int reorder, MPI_Comm* comm_cart);
int MPI_Cart_get(MPI_Comm comm, int maxdims, int* dims,
        int* periods, int* coords);
int MPI_Cart_map(MPI_Comm comm, int ndims, int* dims,
        int* periods, int* newrank);
int MPI_Cart_rank(MPI_Comm comm, int* coords, int* rank);
int MPI_Cart_shift(MPI_Comm comm, int direction, int disp,
        int* rank_source, int* rank_dest);
int MPI_Cart_sub(MPI_Comm comm, int* remain_dims, MPI_Comm* new_comm);
int MPI_Cartdim_get(MPI_Comm comm, int* ndims);
int MPI_Close_port(char* port_name);
int MPI_Comm_accept(char* port_name, MPI_Info info, int root,
        MPI_Comm comm, MPI_Comm* newcomm);
MPI_Fint MPI_Comm_c2f(MPI_Comm comm);
int MPI_Comm_call_errhandler(MPI_Comm comm, int errorcode);
int MPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int* result);
int MPI_Comm_connect(char* port_name, MPI_Info info, int root,
        MPI_Comm comm, MPI_Comm* newcomm);
int MPI_Comm_create_errhandler(MPI_Comm_errhandler_function* function_,
        MPI_Errhandler* errhandler);
int MPI_Comm_create_keyval(MPI_Comm_copy_attr_function* comm_copy_attr_fn,
        MPI_Comm_delete_attr_function* comm_delete_attr_fn,
        int* comm_keyval, void* extra_state);
int MPI_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm* newcomm);
int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm* newcomm);
int MPI_Comm_delete_attr(MPI_Comm comm, int comm_keyval);
int MPI_Comm_disconnect(MPI_Comm* comm);
int MPI_Comm_dup(MPI_Comm comm, MPI_Comm* newcomm);
int MPI_Comm_idup(MPI_Comm comm, MPI_Comm* newcomm, MPI_Request* request);
int MPI_Comm_dup_with_info(MPI_Comm comm, MPI_Info info, MPI_Comm* newcomm);
MPI_Comm MPI_Comm_f2c(MPI_Fint comm);
int MPI_Comm_free_keyval(int* comm_keyval);
int MPI_Comm_free(MPI_Comm* comm);
int MPI_Comm_get_attr(MPI_Comm comm, int comm_keyval,
        void* attribute_val, int* flag);
int MPI_Dist_graph_create(MPI_Comm comm_old, int n, int* nodes,
        int* degrees, int* targets,
        int* weights, MPI_Info info,
        int reorder, MPI_Comm*  newcomm);
int MPI_Dist_graph_create_adjacent(MPI_Comm comm_old,
        int indegree, int* sources,
        int* sourceweights,
        int outdegree,
        int* destinations,
        int* destweights,
        MPI_Info info, int reorder,
        MPI_Comm* comm_dist_graph);
int MPI_Dist_graph_neighbors(MPI_Comm comm, int maxindegree,
        int* sources, int* sourceweights,
        int maxoutdegree,
        int* destinations,
        int* destweights);
int MPI_Dist_graph_neighbors_count(MPI_Comm comm,
        int* inneighbors,
        int* outneighbors,
        int* weighted);
int MPI_Comm_get_errhandler(MPI_Comm comm, MPI_Errhandler* erhandler);
int MPI_Comm_get_info(MPI_Comm comm, MPI_Info* info_used);
int MPI_Comm_get_name(MPI_Comm comm, char* comm_name, int* resultlen);
int MPI_Comm_get_parent(MPI_Comm* parent);
int MPI_Comm_group(MPI_Comm comm, MPI_Group* group);
int MPI_Comm_join(int fd, MPI_Comm* intercomm);
int MPI_Comm_rank(MPI_Comm comm, int* rank);
int MPI_Comm_remote_group(MPI_Comm comm, MPI_Group* group);
int MPI_Comm_remote_size(MPI_Comm comm, int* size);
int MPI_Comm_set_attr(MPI_Comm comm, int comm_keyval, void* attribute_val);
int MPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler);
int MPI_Comm_set_info(MPI_Comm comm, MPI_Info info);
int MPI_Comm_set_name(MPI_Comm comm, char* comm_name);
int MPI_Comm_size(MPI_Comm comm, int* size);
int MPI_Comm_spawn(char* command, char** argv, int maxprocs, MPI_Info info,
        int root, MPI_Comm comm, MPI_Comm* intercomm,
        int* array_of_errcodes);
int MPI_Comm_spawn_multiple(int count, char** array_of_commands, char*** array_of_argv,
        int* array_of_maxprocs, MPI_Info* array_of_info,
        int root, MPI_Comm comm, MPI_Comm* intercomm,
        int* array_of_errcodes);
int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm* newcomm);
int MPI_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm* newcomm);
int MPI_Comm_test_inter(MPI_Comm comm, int* flag);
int MPI_Compare_and_swap(void* origin_addr, void* compare_addr,
        void* result_addr, MPI_Datatype datatype, int target_rank,
        MPI_Aint target_disp, MPI_Win win);
int MPI_Dims_create(int nnodes, int ndims, int* dims);
MPI_Fint MPI_Errhandler_c2f(MPI_Errhandler errhandler);
MPI_Errhandler MPI_Errhandler_f2c(MPI_Fint errhandler);
int MPI_Errhandler_free(MPI_Errhandler* errhandler);
int MPI_Error_class(int errorcode, int* errorclass);
int MPI_Error_string(int errorcode, char* string, int* resultlen);
int MPI_Exscan(void* sendbuf, void* recvbuf, int count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPI_Fetch_and_op(void* origin_addr, void* result_addr, MPI_Datatype datatype,
        int target_rank, MPI_Aint target_disp, MPI_Op op, MPI_Win win);
int MPI_Iexscan(void* sendbuf, void* recvbuf, int count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request* request);

MPI_Fint MPI_File_c2f(MPI_File file);
MPI_File MPI_File_f2c(MPI_Fint file);
int MPI_File_call_errhandler(MPI_File fh, int errorcode);
int MPI_File_create_errhandler(MPI_File_errhandler_function* function_,
                               MPI_Errhandler* errhandler);
int MPI_File_set_errhandler( MPI_File file, MPI_Errhandler errhandler);
int MPI_File_get_errhandler( MPI_File file, MPI_Errhandler* errhandler);
int MPI_File_open(MPI_Comm comm, char* filename, int amode,
                  MPI_Info info, MPI_File* fh);
int MPI_File_close(MPI_File* fh);
int MPI_File_delete(char* filename, MPI_Info info);
int MPI_File_set_size(MPI_File fh, MPI_Offset size);
int MPI_File_preallocate(MPI_File fh, MPI_Offset size);
int MPI_File_get_size(MPI_File fh, MPI_Offset* size);
int MPI_File_get_group(MPI_File fh, MPI_Group* group);
int MPI_File_get_amode(MPI_File fh, int* amode);
int MPI_File_set_info(MPI_File fh, MPI_Info info);
int MPI_File_get_info(MPI_File fh, MPI_Info* info_used);
int MPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype,
                      MPI_Datatype filetype, char* datarep, MPI_Info info);
int MPI_File_get_view(MPI_File fh, MPI_Offset* disp,
                      MPI_Datatype* etype,
                      MPI_Datatype* filetype, char* datarep);
int MPI_File_read_at(MPI_File fh, MPI_Offset offset, void* buf,
                     int count, MPI_Datatype datatype, MPI_Status* status);
int MPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void* buf,
                         int count, MPI_Datatype datatype, MPI_Status* status);
int MPI_File_write_at(MPI_File fh, MPI_Offset offset, void* buf,
                      int count, MPI_Datatype datatype, MPI_Status* status);
int MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, void* buf,
                          int count, MPI_Datatype datatype, MPI_Status* status);
int MPI_File_iread_at(MPI_File fh, MPI_Offset offset, void* buf,
                      int count, MPI_Datatype datatype, MPI_Request* request);
int MPI_File_iwrite_at(MPI_File fh, MPI_Offset offset, void* buf,
                       int count, MPI_Datatype datatype, MPI_Request* request);
int MPI_File_read(MPI_File fh, void* buf, int count,
                  MPI_Datatype datatype, MPI_Status* status);
int MPI_File_read_all(MPI_File fh, void* buf, int count,
                      MPI_Datatype datatype, MPI_Status* status);
int MPI_File_write(MPI_File fh, void* buf, int count,
                   MPI_Datatype datatype, MPI_Status* status);
int MPI_File_write_all(MPI_File fh, void* buf, int count,
                       MPI_Datatype datatype, MPI_Status* status);
int MPI_File_iread(MPI_File fh, void* buf, int count,
                   MPI_Datatype datatype, MPI_Request* request);
int MPI_File_iwrite(MPI_File fh, void* buf, int count,
                    MPI_Datatype datatype, MPI_Request* request);
int MPI_File_seek(MPI_File fh, MPI_Offset offset, int whence);
int MPI_File_get_position(MPI_File fh, MPI_Offset* offset);
int MPI_File_get_byte_offset(MPI_File fh, MPI_Offset offset,
                             MPI_Offset* disp);
int MPI_File_read_shared(MPI_File fh, void* buf, int count,
                         MPI_Datatype datatype, MPI_Status* status);
int MPI_File_write_shared(MPI_File fh, void* buf, int count,
                          MPI_Datatype datatype, MPI_Status* status);
int MPI_File_iread_shared(MPI_File fh, void* buf, int count,
                          MPI_Datatype datatype, MPI_Request* request);
int MPI_File_iwrite_shared(MPI_File fh, void* buf, int count,
                           MPI_Datatype datatype, MPI_Request* request);
int MPI_File_read_ordered(MPI_File fh, void* buf, int count,
                          MPI_Datatype datatype, MPI_Status* status);
int MPI_File_write_ordered(MPI_File fh, void* buf, int count,
                           MPI_Datatype datatype, MPI_Status* status);
int MPI_File_seek_shared(MPI_File fh, MPI_Offset offset, int whence);
int MPI_File_get_position_shared(MPI_File fh, MPI_Offset* offset);
int MPI_File_read_at_all_begin(MPI_File fh, MPI_Offset offset, void* buf,
                               int count, MPI_Datatype datatype);
int MPI_File_read_at_all_end(MPI_File fh, void* buf, MPI_Status* status);
int MPI_File_write_at_all_begin(MPI_File fh, MPI_Offset offset, void* buf,
                                int count, MPI_Datatype datatype);
int MPI_File_write_at_all_end(MPI_File fh, void* buf, MPI_Status* status);
int MPI_File_read_all_begin(MPI_File fh, void* buf, int count,
                            MPI_Datatype datatype);
int MPI_File_read_all_end(MPI_File fh, void* buf, MPI_Status* status);
int MPI_File_write_all_begin(MPI_File fh, void* buf, int count,
                             MPI_Datatype datatype);
int MPI_File_write_all_end(MPI_File fh, void* buf, MPI_Status* status);
int MPI_File_read_ordered_begin(MPI_File fh, void* buf, int count,
                                MPI_Datatype datatype);
int MPI_File_read_ordered_end(MPI_File fh, void* buf, MPI_Status* status);
int MPI_File_write_ordered_begin(MPI_File fh, void* buf, int count,
                                 MPI_Datatype datatype);
int MPI_File_write_ordered_end(MPI_File fh, void* buf, MPI_Status* status);
int MPI_File_get_type_extent(MPI_File fh, MPI_Datatype datatype,
                             MPI_Aint* extent);
int MPI_File_set_atomicity(MPI_File fh, int flag);
int MPI_File_get_atomicity(MPI_File fh, int* flag);
int MPI_File_sync(MPI_File fh);

int MPI_Finalize();
int MPI_Finalized(int* flag);
int MPI_Free_mem(void* base);
int MPI_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount, MPI_Datatype recvtype,
        int root, MPI_Comm comm);
int MPI_Igather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount, MPI_Datatype recvtype,
        int root, MPI_Comm comm, MPI_Request* request);
int MPI_Gatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int* recvcounts, int* displs,
        MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPI_Igatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int* recvcounts, int* displs,
        MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request* request);
int MPI_Get_address(void* location, MPI_Aint* address);
int MPI_Get_count(MPI_Status* status, MPI_Datatype datatype, int* count);
int MPI_Get_elements(MPI_Status* status, MPI_Datatype datatype, int* count);
static if(OMPI_MAJOR_VERSION == 1 && OMPI_MINOR_VERSION >= 8)
    int MPI_Get_elements_x(MPI_Status* status, MPI_Datatype datatype, MPI_Count* count);
int MPI_Get(void* origin_addr, int origin_count,
        MPI_Datatype origin_datatype, int target_rank,
        MPI_Aint target_disp, int target_count,
        MPI_Datatype target_datatype, MPI_Win win);
int MPI_Get_accumulate(void* origin_addr, int origin_count, MPI_Datatype origin_datatype,
        void* result_addr, int result_count, MPI_Datatype result_datatype,
        int target_rank, MPI_Aint target_disp, int target_count,
        MPI_Datatype target_datatype, MPI_Op op, MPI_Win win);
int MPI_Get_library_version(char* version_, int* resultlen);
int MPI_Get_processor_name(char* name, int* resultlen);
int MPI_Get_version(int* version_, int* subversion);
int MPI_Graph_create(MPI_Comm comm_old, int nnodes, int* index,
        int* edges, int reorder, MPI_Comm* comm_graph);
int MPI_Graph_get(MPI_Comm comm, int maxindex, int maxedges,
        int* index, int* edges);
int MPI_Graph_map(MPI_Comm comm, int nnodes, int* index, int* edges,
        int* newrank);
int MPI_Graph_neighbors_count(MPI_Comm comm, int rank, int* nneighbors);
int MPI_Graph_neighbors(MPI_Comm comm, int rank, int maxneighbors,
        int* neighbors);
int MPI_Graphdims_get(MPI_Comm comm, int* nnodes, int* nedges);
int MPI_Grequest_complete(MPI_Request request);
int MPI_Grequest_start(MPI_Grequest_query_function* query_fn,
        MPI_Grequest_free_function* free_fn,
        MPI_Grequest_cancel_function* cancel_fn,
        void* extra_state, MPI_Request* request);
MPI_Fint MPI_Group_c2f(MPI_Group group);
int MPI_Group_compare(MPI_Group group1, MPI_Group group2, int* result);
int MPI_Group_difference(MPI_Group group1, MPI_Group group2,
        MPI_Group* newgroup);
int MPI_Group_excl(MPI_Group group, int n, int* ranks,
        MPI_Group* newgroup);
MPI_Group MPI_Group_f2c(MPI_Fint group);
int MPI_Group_free(MPI_Group* group);
int MPI_Group_incl(MPI_Group group, int n, int* ranks,
        MPI_Group* newgroup);
int MPI_Group_intersection(MPI_Group group1, MPI_Group group2,
        MPI_Group* newgroup);
int MPI_Group_range_excl(MPI_Group group, int n, int[3]* ranges,
        MPI_Group* newgroup);
int MPI_Group_range_incl(MPI_Group group, int n, int[3]* ranges,
        MPI_Group* newgroup);
int MPI_Group_rank(MPI_Group group, int* rank);
int MPI_Group_size(MPI_Group group, int* size);
int MPI_Group_translate_ranks(MPI_Group group1, int n, int* ranks1,
        MPI_Group group2, int* ranks2);
int MPI_Group_union(MPI_Group group1, MPI_Group group2,
        MPI_Group* newgroup);
int MPI_Ibsend(void* buf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm, MPI_Request* request);
int MPI_Improbe(int source, int tag, MPI_Comm comm,
                int* flag, MPI_Message* message,
                MPI_Status* status);
int MPI_Imrecv(void* buf, int count, MPI_Datatype type,
               MPI_Message* message, MPI_Request* request);
MPI_Fint MPI_Info_c2f(MPI_Info info);
int MPI_Info_create(MPI_Info* info);
int MPI_Info_delete(MPI_Info info, char* key);
int MPI_Info_dup(MPI_Info info, MPI_Info* newinfo);
MPI_Info MPI_Info_f2c(MPI_Fint info);
int MPI_Info_free(MPI_Info* info);
int MPI_Info_get(MPI_Info info, char* key, int valuelen,
        char* value, int* flag);
int MPI_Info_get_nkeys(MPI_Info info, int* nkeys);
int MPI_Info_get_nthkey(MPI_Info info, int n, char* key);
int MPI_Info_get_valuelen(MPI_Info info, char* key, int* valuelen,
        int* flag);
int MPI_Info_set(MPI_Info info, char* key, char* value);
int MPI_Init(int* argc, char*** argv);
int MPI_Initialized(int* flag);
int MPI_Init_thread(int* argc, char*** argv, int required,
        int* provided);
int MPI_Intercomm_create(MPI_Comm local_comm, int local_leader,
        MPI_Comm bridge_comm, int remote_leader,
        int tag, MPI_Comm* newintercomm);
int MPI_Intercomm_merge(MPI_Comm intercomm, int high,
        MPI_Comm* newintercomm);
int MPI_Iprobe(int source, int tag, MPI_Comm comm, int* flag,
        MPI_Status* status);
int MPI_Irecv(void* buf, int count, MPI_Datatype datatype, int source,
        int tag, MPI_Comm comm, MPI_Request* request);
int MPI_Irsend(void* buf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm, MPI_Request* request);
int MPI_Isend(void* buf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm, MPI_Request* request);
int MPI_Issend(void* buf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm, MPI_Request* request);
int MPI_Is_thread_main(int* flag);
int MPI_Lookup_name(char* service_name, MPI_Info info, char* port_name);
MPI_Fint MPI_Message_c2f(MPI_Message message);
MPI_Message MPI_Message_f2c(MPI_Fint message);
int MPI_Mprobe(int source, int tag, MPI_Comm comm,
        MPI_Message* message,
        MPI_Status* status);
int MPI_Mrecv(void* buf, int count, MPI_Datatype type,
        MPI_Message* message, MPI_Status* status);
int MPI_Neighbor_allgather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount, MPI_Datatype recvtype,
        MPI_Comm comm);
int MPI_Ineighbor_allgather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount, MPI_Datatype recvtype,
        MPI_Comm comm, MPI_Request* request);
int MPI_Neighbor_allgatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int* recvcounts, int* displs,
        MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Ineighbor_allgatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int* recvcounts, int* displs,
        MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request);
int MPI_Neighbor_alltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount, MPI_Datatype recvtype,
        MPI_Comm comm);
int MPI_Ineighbor_alltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount, MPI_Datatype recvtype,
        MPI_Comm comm, MPI_Request* request);
int MPI_Neighbor_alltoallv(void* sendbuf, int* sendcounts, int* sdispls,  MPI_Datatype sendtype,
        void* recvbuf, int* recvcounts, int* rdispls, MPI_Datatype recvtype,
        MPI_Comm comm);
int MPI_Ineighbor_alltoallv(void* sendbuf, int* sendcounts, int* sdispls, MPI_Datatype sendtype,
        void* recvbuf, int* recvcounts, int* rdispls, MPI_Datatype recvtype,
        MPI_Comm comm, MPI_Request* request);
int MPI_Neighbor_alltoallw(void* sendbuf, int* sendcounts, MPI_Aint* sdispls, MPI_Datatype* sendtypes,
        void* recvbuf, int* recvcounts, MPI_Aint* rdispls, MPI_Datatype* recvtypes,
        MPI_Comm comm);
int MPI_Ineighbor_alltoallw(void* sendbuf, int* sendcounts, MPI_Aint* sdispls, MPI_Datatype* sendtypes,
        void* recvbuf, int* recvcounts, MPI_Aint* rdispls, MPI_Datatype* recvtypes,
        MPI_Comm comm, MPI_Request* request);
MPI_Fint MPI_Op_c2f(MPI_Op op);
int MPI_Op_commutative(MPI_Op op, int* commute);
int MPI_Op_create(MPI_User_function* function_, int commute, MPI_Op* op);
int MPI_Open_port(MPI_Info info, char* port_name);
MPI_Op MPI_Op_f2c(MPI_Fint op);
int MPI_Op_free(MPI_Op* op);
int MPI_Pack_external(char* datarep, void* inbuf, int incount,
        MPI_Datatype datatype, void* outbuf,
        MPI_Aint outsize, MPI_Aint* position);
int MPI_Pack_external_size(char* datarep, int incount,
        MPI_Datatype datatype, MPI_Aint* size);
int MPI_Pack(void* inbuf, int incount, MPI_Datatype datatype,
        void* outbuf, int outsize, int* position, MPI_Comm comm);
int MPI_Pack_size(int incount, MPI_Datatype datatype, MPI_Comm comm,
        int* size);
int MPI_Pcontrol(int level, ...);
int MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status* status);
int MPI_Publish_name(char* service_name, MPI_Info info,
        char* port_name);
int MPI_Put(void* origin_addr, int origin_count, MPI_Datatype origin_datatype,
        int target_rank, MPI_Aint target_disp, int target_count,
        MPI_Datatype target_datatype, MPI_Win win);
int MPI_Query_thread(int* provided);
int MPI_Raccumulate(void* origin_addr, int origin_count, MPI_Datatype origin_datatype, 
        int target_rank, MPI_Aint target_disp, int target_count, 
        MPI_Datatype target_datatype, MPI_Op op, MPI_Win win, MPI_Request* request);
int MPI_Recv_init(void* buf, int count, MPI_Datatype datatype, int source,
        int tag, MPI_Comm comm, MPI_Request* request);
int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source,
        int tag, MPI_Comm comm, MPI_Status* status);
int MPI_Reduce(void* sendbuf, void* recvbuf, int count,
        MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
int MPI_Ireduce(void* sendbuf, void* recvbuf, int count,
        MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request* request);
int MPI_Reduce_local(void* inbuf, void* inoutbuf, int count,
        MPI_Datatype datatype, MPI_Op op);
int MPI_Reduce_scatter(void* sendbuf, void* recvbuf, int* recvcounts,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPI_Ireduce_scatter(void* sendbuf, void* recvbuf, int* recvcounts,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request* request);
int MPI_Reduce_scatter_block(void* sendbuf, void* recvbuf, int recvcount,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPI_Ireduce_scatter_block(void* sendbuf, void* recvbuf, int recvcount,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request* request);
int MPI_Register_datarep(char* datarep,
        MPI_Datarep_conversion_function* read_conversion_fn,
        MPI_Datarep_conversion_function* write_conversion_fn,
        MPI_Datarep_extent_function* dtype_file_extent_fn,
        void* extra_state);
MPI_Fint MPI_Request_c2f(MPI_Request request);
MPI_Request MPI_Request_f2c(MPI_Fint request);
int MPI_Request_free(MPI_Request* request);
int MPI_Request_get_status(MPI_Request request, int* flag,
        MPI_Status* status);
int MPI_Rget(void* origin_addr, int origin_count, MPI_Datatype origin_datatype, 
        int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype,
        MPI_Win win, MPI_Request* request);
int MPI_Rget_accumulate(void* origin_addr, int origin_count, MPI_Datatype origin_datatype,
        void* result_addr, int result_count, MPI_Datatype result_datatype,
        int target_rank, MPI_Aint target_disp, int target_count, 
        MPI_Datatype target_datatype, MPI_Op op,
        MPI_Win win, MPI_Request* request);
int MPI_Rput(void* origin_addr, int origin_count, MPI_Datatype origin_datatype,
        int target_rank, MPI_Aint target_disp, int target_cout, 
        MPI_Datatype target_datatype, MPI_Win win, MPI_Request* request);
int MPI_Rsend(void* ibuf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm);
int MPI_Rsend_init(void* buf, int count, MPI_Datatype datatype,
        int dest, int tag, MPI_Comm comm,
        MPI_Request* request);
int MPI_Scan(void* sendbuf, void* recvbuf, int count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPI_Iscan(void* sendbuf, void* recvbuf, int count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request* request);
int MPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount, MPI_Datatype recvtype,
        int root, MPI_Comm comm);
int MPI_Iscatter(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount, MPI_Datatype recvtype,
        int root, MPI_Comm comm, MPI_Request* request);
int MPI_Scatterv(void* sendbuf, int* sendcounts, int* displs,
        MPI_Datatype sendtype, void* recvbuf, int recvcount,
        MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPI_Iscatterv(void* sendbuf, int* sendcounts, int* displs,
        MPI_Datatype sendtype, void* recvbuf, int recvcount,
        MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request* request);
int MPI_Send_init(void* buf, int count, MPI_Datatype datatype,
        int dest, int tag, MPI_Comm comm,
        MPI_Request* request);
int MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm);
int MPI_Sendrecv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        int dest, int sendtag, void* recvbuf, int recvcount,
        MPI_Datatype recvtype, int source, int recvtag,
        MPI_Comm comm,  MPI_Status* status);
int MPI_Sendrecv_replace(void*  buf, int count, MPI_Datatype datatype,
        int dest, int sendtag, int source, int recvtag,
        MPI_Comm comm, MPI_Status* status);
int MPI_Ssend_init(void* buf, int count, MPI_Datatype datatype,
        int dest, int tag, MPI_Comm comm,
        MPI_Request* request);
int MPI_Ssend(void* buf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm);
int MPI_Start(MPI_Request* request);
int MPI_Startall(int count, MPI_Request* array_of_requests);
int MPI_Status_c2f(MPI_Status* c_status, MPI_Fint* f_status);
int MPI_Status_f2c(MPI_Fint* f_status, MPI_Status* c_status);
int MPI_Status_set_cancelled(MPI_Status* status, int flag);
int MPI_Status_set_elements(MPI_Status* status, MPI_Datatype datatype,
        int count);
int MPI_Status_set_elements_x(MPI_Status* status, MPI_Datatype datatype,
        MPI_Count count);
int MPI_Testall(int count, MPI_Request* array_of_requests, int* flag,
        MPI_Status* array_of_statuses);
int MPI_Testany(int count, MPI_Request* array_of_requests, int* index,
        int* flag, MPI_Status* status);
int MPI_Test(MPI_Request* request, int* flag, MPI_Status* status);
int MPI_Test_cancelled(MPI_Status* status, int* flag);
int MPI_Testsome(int incount, MPI_Request* array_of_requests,
        int* outcount, int* array_of_indices,
        MPI_Status* array_of_statuses);
int MPI_Topo_test(MPI_Comm comm, int* status);
MPI_Fint MPI_Type_c2f(MPI_Datatype datatype);
int MPI_Type_commit(MPI_Datatype* type);
int MPI_Type_contiguous(int count, MPI_Datatype oldtype,
        MPI_Datatype* newtype);
int MPI_Type_create_darray(int size, int rank, int ndims,
        int* gsize_array, int* distrib_array,
        int* darg_array, int* psize_array,
        int order, MPI_Datatype oldtype,
        MPI_Datatype* newtype);
int MPI_Type_create_f90_complex(int p, int r, MPI_Datatype* newtype);
int MPI_Type_create_f90_integer(int r, MPI_Datatype* newtype);
int MPI_Type_create_f90_real(int p, int r, MPI_Datatype* newtype);
int MPI_Type_create_hindexed_block(int count, int blocklength,
        MPI_Aint* array_of_displacements,
        MPI_Datatype oldtype,
        MPI_Datatype* newtype);
int MPI_Type_create_hindexed(int count, int* array_of_blocklengths,
        MPI_Aint* array_of_displacements,
        MPI_Datatype oldtype,
        MPI_Datatype* newtype);
int MPI_Type_create_hvector(int count, int blocklength, MPI_Aint stride,
        MPI_Datatype oldtype,
        MPI_Datatype* newtype);
int MPI_Type_create_keyval(MPI_Type_copy_attr_function* type_copy_attr_fn,
        MPI_Type_delete_attr_function* type_delete_attr_fn,
        int* type_keyval, void* extra_state);
int MPI_Type_create_indexed_block(int count, int blocklength,
        int* array_of_displacements,
        MPI_Datatype oldtype,
        MPI_Datatype* newtype);
int MPI_Type_create_struct(int count, int* array_of_block_lengths,
        MPI_Aint* array_of_displacements,
        MPI_Datatype* array_of_types,
        MPI_Datatype* newtype);
int MPI_Type_create_subarray(int ndims, int* size_array, int* subsize_array,
        int* start_array, int order,
        MPI_Datatype oldtype, MPI_Datatype* newtype);
int MPI_Type_create_resized(MPI_Datatype oldtype, MPI_Aint lb,
        MPI_Aint extent, MPI_Datatype* newtype);
int MPI_Type_delete_attr(MPI_Datatype type, int type_keyval);
int MPI_Type_dup(MPI_Datatype type, MPI_Datatype* newtype);
int MPI_Type_free(MPI_Datatype* type);
int MPI_Type_free_keyval(int* type_keyval);
MPI_Datatype MPI_Type_f2c(MPI_Fint datatype);
int MPI_Type_get_attr(MPI_Datatype type, int type_keyval,
        void* attribute_val, int* flag);
int MPI_Type_get_contents(MPI_Datatype mtype, int max_integers,
        int max_addresses, int max_datatypes,
        int* array_of_integers,
        MPI_Aint* array_of_addresses,
        MPI_Datatype* array_of_datatypes);
int MPI_Type_get_envelope(MPI_Datatype type, int* num_integers,
        int* num_addresses, int* num_datatypes,
        int* combiner);
int MPI_Type_get_extent(MPI_Datatype type, MPI_Aint* lb,
        MPI_Aint* extent);
int MPI_Type_get_extent_x(MPI_Datatype type, MPI_Count* lb,
        MPI_Count* extent);
int MPI_Type_get_name(MPI_Datatype type, char* type_name,
        int* resultlen);
int MPI_Type_get_true_extent(MPI_Datatype datatype, MPI_Aint* true_lb,
        MPI_Aint* true_extent);
int MPI_Type_get_true_extent_x(MPI_Datatype datatype, MPI_Count* true_lb,
        MPI_Count* true_extent);
int MPI_Type_indexed(int count, int* array_of_blocklengths,
        int* array_of_displacements,
        MPI_Datatype oldtype, MPI_Datatype* newtype);
int MPI_Type_match_size(int typeclass, int size, MPI_Datatype* type);
int MPI_Type_set_attr(MPI_Datatype type, int type_keyval,
        void* attr_val);
int MPI_Type_set_name(MPI_Datatype type, char* type_name);
int MPI_Type_size(MPI_Datatype type, int* size);
int MPI_Type_size_x(MPI_Datatype type, MPI_Count* size);
int MPI_Type_vector(int count, int blocklength, int stride,
        MPI_Datatype oldtype, MPI_Datatype* newtype);
int MPI_Unpack(void* inbuf, int insize, int* position,
        void* outbuf, int outcount, MPI_Datatype datatype,
        MPI_Comm comm);
int MPI_Unpublish_name(char* service_name, MPI_Info info, char* port_name);
int MPI_Unpack_external (char* datarep, void* inbuf, MPI_Aint insize,
        MPI_Aint* position, void* outbuf, int outcount,
        MPI_Datatype datatype);
int MPI_Waitall(int count, MPI_Request* array_of_requests,
        MPI_Status* array_of_statuses);
int MPI_Waitany(int count, MPI_Request* array_of_requests,
        int* index, MPI_Status* status);
int MPI_Wait(MPI_Request* request, MPI_Status* status);
int MPI_Waitsome(int incount, MPI_Request* array_of_requests,
        int* outcount, int* array_of_indices,
        MPI_Status* array_of_statuses);
int MPI_Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info,
        MPI_Comm comm, void* baseptr, MPI_Win* win);
int MPI_Win_allocate_shared(MPI_Aint size, int disp_unit, MPI_Info info,
        MPI_Comm comm, void* baseptr, MPI_Win* win);
int MPI_Win_attach(MPI_Win win, void* base, MPI_Aint size);
MPI_Fint MPI_Win_c2f(MPI_Win win);
int MPI_Win_call_errhandler(MPI_Win win, int errorcode);
int MPI_Win_complete(MPI_Win win);
int MPI_Win_create(void* base, MPI_Aint size, int disp_unit,
        MPI_Info info, MPI_Comm comm, MPI_Win* win);
int MPI_Win_create_dynamic(MPI_Info info, MPI_Comm comm, MPI_Win* win);
int MPI_Win_create_errhandler(MPI_Win_errhandler_function* function_,
        MPI_Errhandler* errhandler);
int MPI_Win_create_keyval(MPI_Win_copy_attr_function* win_copy_attr_fn,
        MPI_Win_delete_attr_function* win_delete_attr_fn,
        int* win_keyval, void* extra_state);
int MPI_Win_delete_attr(MPI_Win win, int win_keyval);
int MPI_Win_detach(MPI_Win win, void* base);
MPI_Win MPI_Win_f2c(MPI_Fint win);
int MPI_Win_fence(int assert_, MPI_Win win);
int MPI_Win_flush(int rank, MPI_Win win);
int MPI_Win_flush_all(MPI_Win win);
int MPI_Win_flush_local(int rank, MPI_Win win);
int MPI_Win_flush_local_all(MPI_Win win);
int MPI_Win_free(MPI_Win* win);
int MPI_Win_free_keyval(int* win_keyval);
int MPI_Win_get_attr(MPI_Win win, int win_keyval,
        void* attribute_val, int* flag);
int MPI_Win_get_errhandler(MPI_Win win, MPI_Errhandler* errhandler);
int MPI_Win_get_group(MPI_Win win, MPI_Group* group);
int MPI_Win_get_info(MPI_Win win, MPI_Info* info_used);
int MPI_Win_get_name(MPI_Win win, char* win_name, int* resultlen);
int MPI_Win_lock(int lock_type, int rank, int assert_, MPI_Win win);
int MPI_Win_lock_all(int assert_, MPI_Win win);
int MPI_Win_post(MPI_Group group, int assert_, MPI_Win win);
int MPI_Win_set_attr(MPI_Win win, int win_keyval, void* attribute_val);
int MPI_Win_set_errhandler(MPI_Win win, MPI_Errhandler errhandler);
int MPI_Win_set_info(MPI_Win win, MPI_Info info);
int MPI_Win_set_name(MPI_Win win, char* win_name);
int MPI_Win_shared_query(MPI_Win win, int rank, MPI_Aint* size, int* disp_unit, void* baseptr);
int MPI_Win_start(MPI_Group group, int assert_, MPI_Win win);
int MPI_Win_sync(MPI_Win win);
int MPI_Win_test(MPI_Win win, int* flag);
int MPI_Win_unlock(int rank, MPI_Win win);
int MPI_Win_unlock_all(MPI_Win win);
int MPI_Win_wait(MPI_Win win);
double MPI_Wtick();
double MPI_Wtime();

/*
 * Profiling MPI API
 */
int PMPI_Abort(MPI_Comm comm, int errorcode);
int PMPI_Accumulate(void* origin_addr, int origin_count, MPI_Datatype origin_datatype,
        int target_rank, MPI_Aint target_disp, int target_count,
        MPI_Datatype target_datatype, MPI_Op op, MPI_Win win);
int PMPI_Add_error_class(int* errorclass);
int PMPI_Add_error_code(int errorclass, int* errorcode);
int PMPI_Add_error_string(int errorcode, char* string);
int PMPI_Allgather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount,
        MPI_Datatype recvtype, MPI_Comm comm);
int PMPI_Iallgather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount,
        MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request);
int PMPI_Allgatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int* recvcounts,
        int* displs, MPI_Datatype recvtype, MPI_Comm comm);
int PMPI_Iallgatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int* recvcounts,
        int* displs, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request);
int PMPI_Alloc_mem(MPI_Aint size, MPI_Info info,
        void* baseptr);
int PMPI_Allreduce(void* sendbuf, void* recvbuf, int count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int PMPI_Iallreduce(void* sendbuf, void* recvbuf, int count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request* request);
int PMPI_Alltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount,
        MPI_Datatype recvtype, MPI_Comm comm);
int PMPI_Ialltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount,
        MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request);
int PMPI_Alltoallv(void* sendbuf, int* sendcounts, int* sdispls,
        MPI_Datatype sendtype, void* recvbuf, int* recvcounts,
        int* rdispls, MPI_Datatype recvtype, MPI_Comm comm);
int PMPI_Ialltoallv(void* sendbuf, int* sendcounts, int* sdispls,
        MPI_Datatype sendtype, void* recvbuf, int* recvcounts,
        int* rdispls, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request);
int PMPI_Alltoallw(void* sendbuf, int* sendcounts, int* sdispls, MPI_Datatype* sendtypes,
        void* recvbuf, int* recvcounts, int* rdispls, MPI_Datatype* recvtypes,
        MPI_Comm comm);
int PMPI_Ialltoallw(void* sendbuf, int* sendcounts, int* sdispls, MPI_Datatype* sendtypes,
        void* recvbuf, int* recvcounts, int* rdispls, MPI_Datatype* recvtypes,
        MPI_Comm comm, MPI_Request* request);
int PMPI_Dist_graph_create(MPI_Comm comm_old, int n, int* nodes,
        int* degrees, int* targets,
        int* weights, MPI_Info info,
        int reorder, MPI_Comm*  newcomm);
int PMPI_Dist_graph_create_adjacent(MPI_Comm comm_old,
        int indegree, int* sources,
        int* sourceweights,
        int outdegree,
        int* destinations,
        int* destweights,
        MPI_Info info, int reorder,
        MPI_Comm* comm_dist_graph);
int PMPI_Dist_graph_neighbors(MPI_Comm comm, int maxindegree,
        int* sources, int* sourceweights,
        int maxoutdegree,
        int* destinations,
        int* destweights);
int PMPI_Dist_graph_neighbors_count(MPI_Comm comm,
        int* inneighbors,
        int* outneighbors,
        int* weighted);
int PMPI_Barrier(MPI_Comm comm);
int PMPI_Ibarrier(MPI_Comm comm, MPI_Request* request);
int PMPI_Bcast(void* buffer, int count, MPI_Datatype datatype,
        int root, MPI_Comm comm);
int PMPI_Ibcast(void* buffer, int count, MPI_Datatype datatype,
        int root, MPI_Comm comm,
        MPI_Request* request);
int PMPI_Bsend(void* buf, int count, MPI_Datatype datatype,
        int dest, int tag, MPI_Comm comm);
int PMPI_Bsend_init(void* buf, int count, MPI_Datatype datatype,
        int dest, int tag, MPI_Comm comm, MPI_Request* request);
int PMPI_Buffer_attach(void* buffer, int size);
int PMPI_Buffer_detach(void* buffer, int* size);
int PMPI_Cancel(MPI_Request* request);
int PMPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int* coords);
int PMPI_Cart_create(MPI_Comm old_comm, int ndims, int* dims,
        int* periods, int reorder, MPI_Comm* comm_cart);
int PMPI_Cart_get(MPI_Comm comm, int maxdims, int* dims,
        int* periods, int* coords);
int PMPI_Cart_map(MPI_Comm comm, int ndims, int* dims,
        int* periods, int* newrank);
int PMPI_Cart_rank(MPI_Comm comm, int* coords, int* rank);
int PMPI_Cart_shift(MPI_Comm comm, int direction, int disp,
        int* rank_source, int* rank_dest);
int PMPI_Cart_sub(MPI_Comm comm, int* remain_dims, MPI_Comm* new_comm);
int PMPI_Cartdim_get(MPI_Comm comm, int* ndims);
int PMPI_Close_port(char* port_name);
int PMPI_Comm_accept(char* port_name, MPI_Info info, int root,
        MPI_Comm comm, MPI_Comm* newcomm);
MPI_Fint PMPI_Comm_c2f(MPI_Comm comm);
int PMPI_Comm_call_errhandler(MPI_Comm comm, int errorcode);
int PMPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int* result);
int PMPI_Comm_connect(char* port_name, MPI_Info info, int root,
        MPI_Comm comm, MPI_Comm* newcomm);
int PMPI_Comm_create_errhandler(MPI_Comm_errhandler_function* function_,
        MPI_Errhandler* errhandler);
int PMPI_Comm_create_keyval(MPI_Comm_copy_attr_function* comm_copy_attr_fn,
        MPI_Comm_delete_attr_function* comm_delete_attr_fn,
        int* comm_keyval, void* extra_state);
int PMPI_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm* newcomm);
int PMPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm* newcomm);
int PMPI_Comm_delete_attr(MPI_Comm comm, int comm_keyval);
int PMPI_Comm_disconnect(MPI_Comm* comm);
int PMPI_Comm_dup(MPI_Comm comm, MPI_Comm* newcomm);
int PMPI_Comm_idup(MPI_Comm comm, MPI_Comm* newcomm, MPI_Request* request);
int PMPI_Comm_dup_with_info(MPI_Comm comm, MPI_Info info, MPI_Comm* newcomm);
MPI_Comm PMPI_Comm_f2c(MPI_Fint comm);
int PMPI_Comm_free_keyval(int* comm_keyval);
int PMPI_Comm_free(MPI_Comm* comm);
int PMPI_Comm_get_attr(MPI_Comm comm, int comm_keyval,
        void* attribute_val, int* flag);
int PMPI_Comm_get_errhandler(MPI_Comm comm, MPI_Errhandler* erhandler);
int PMPI_Comm_get_info(MPI_Comm comm, MPI_Info* info_used);
int PMPI_Comm_get_name(MPI_Comm comm, char* comm_name, int* resultlen);
int PMPI_Comm_get_parent(MPI_Comm* parent);
int PMPI_Comm_group(MPI_Comm comm, MPI_Group* group);
int PMPI_Comm_join(int fd, MPI_Comm* intercomm);
int PMPI_Comm_rank(MPI_Comm comm, int* rank);
int PMPI_Comm_remote_group(MPI_Comm comm, MPI_Group* group);
int PMPI_Comm_remote_size(MPI_Comm comm, int* size);
int PMPI_Comm_set_attr(MPI_Comm comm, int comm_keyval, void* attribute_val);
int PMPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler);
int PMPI_Comm_set_info(MPI_Comm comm, MPI_Info info);
int PMPI_Comm_set_name(MPI_Comm comm, char* comm_name);
int PMPI_Comm_size(MPI_Comm comm, int* size);
int PMPI_Comm_spawn(char* command, char** argv, int maxprocs, MPI_Info info,
        int root, MPI_Comm comm, MPI_Comm* intercomm,
        int* array_of_errcodes);
int PMPI_Comm_spawn_multiple(int count, char** array_of_commands, char*** array_of_argv,
        int* array_of_maxprocs, MPI_Info* array_of_info,
        int root, MPI_Comm comm, MPI_Comm* intercomm,
        int* array_of_errcodes);
int PMPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm* newcomm);
int PMPI_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm* newcomm);
int PMPI_Comm_test_inter(MPI_Comm comm, int* flag);
int PMPI_Compare_and_swap(void* origin_addr, void* compare_addr,
        void* result_addr, MPI_Datatype datatype, int target_rank,
        MPI_Aint target_disp, MPI_Win win);
int PMPI_Dims_create(int nnodes, int ndims, int* dims);
MPI_Fint PMPI_Errhandler_c2f(MPI_Errhandler errhandler);
MPI_Errhandler PMPI_Errhandler_f2c(MPI_Fint errhandler);
int PMPI_Errhandler_free(MPI_Errhandler* errhandler);
int PMPI_Error_class(int errorcode, int* errorclass);
int PMPI_Error_string(int errorcode, char* string, int* resultlen);
int PMPI_Exscan(void* sendbuf, void* recvbuf, int count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int PMPI_Fetch_and_op(void* origin_addr, void* result_addr, MPI_Datatype datatype,
        int target_rank, MPI_Aint target_disp, MPI_Op op, MPI_Win win);
int PMPI_Iexscan(void* sendbuf, void* recvbuf, int count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request* request);
MPI_Fint PMPI_File_c2f(MPI_File file);
MPI_File PMPI_File_f2c(MPI_Fint file);
int PMPI_File_call_errhandler(MPI_File fh, int errorcode);
int PMPI_File_create_errhandler(MPI_File_errhandler_function* function_,
        MPI_Errhandler* errhandler);
int PMPI_File_set_errhandler( MPI_File file, MPI_Errhandler errhandler);
int PMPI_File_get_errhandler( MPI_File file, MPI_Errhandler* errhandler);
int PMPI_File_open(MPI_Comm comm, char* filename, int amode,
        MPI_Info info, MPI_File* fh);
int PMPI_File_close(MPI_File* fh);
int PMPI_File_delete(char* filename, MPI_Info info);
int PMPI_File_set_size(MPI_File fh, MPI_Offset size);
int PMPI_File_preallocate(MPI_File fh, MPI_Offset size);
int PMPI_File_get_size(MPI_File fh, MPI_Offset* size);
int PMPI_File_get_group(MPI_File fh, MPI_Group* group);
int PMPI_File_get_amode(MPI_File fh, int* amode);
int PMPI_File_set_info(MPI_File fh, MPI_Info info);
int PMPI_File_get_info(MPI_File fh, MPI_Info* info_used);
int PMPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype,
        MPI_Datatype filetype, char* datarep, MPI_Info info);
int PMPI_File_get_view(MPI_File fh, MPI_Offset* disp,
        MPI_Datatype* etype,
        MPI_Datatype* filetype, char* datarep);
int PMPI_File_read_at(MPI_File fh, MPI_Offset offset, void* buf,
        int count, MPI_Datatype datatype, MPI_Status* status);
int PMPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void* buf,
        int count, MPI_Datatype datatype, MPI_Status* status);
int PMPI_File_write_at(MPI_File fh, MPI_Offset offset, void* buf,
        int count, MPI_Datatype datatype, MPI_Status* status);
int PMPI_File_write_at_all(MPI_File fh, MPI_Offset offset, void* buf,
        int count, MPI_Datatype datatype, MPI_Status* status);
int PMPI_File_iread_at(MPI_File fh, MPI_Offset offset, void* buf,
        int count, MPI_Datatype datatype, MPI_Request* request);
int PMPI_File_iwrite_at(MPI_File fh, MPI_Offset offset, void* buf,
        int count, MPI_Datatype datatype, MPI_Request* request);
int PMPI_File_read(MPI_File fh, void* buf, int count,
        MPI_Datatype datatype, MPI_Status* status);
int PMPI_File_read_all(MPI_File fh, void* buf, int count,
        MPI_Datatype datatype, MPI_Status* status);
int PMPI_File_write(MPI_File fh, void* buf, int count,
        MPI_Datatype datatype, MPI_Status* status);
int PMPI_File_write_all(MPI_File fh, void* buf, int count,
        MPI_Datatype datatype, MPI_Status* status);
int PMPI_File_iread(MPI_File fh, void* buf, int count,
        MPI_Datatype datatype, MPI_Request* request);
int PMPI_File_iwrite(MPI_File fh, void* buf, int count,
        MPI_Datatype datatype, MPI_Request* request);
int PMPI_File_seek(MPI_File fh, MPI_Offset offset, int whence);
int PMPI_File_get_position(MPI_File fh, MPI_Offset* offset);
int PMPI_File_get_byte_offset(MPI_File fh, MPI_Offset offset,
        MPI_Offset* disp);
int PMPI_File_read_shared(MPI_File fh, void* buf, int count,
        MPI_Datatype datatype, MPI_Status* status);
int PMPI_File_write_shared(MPI_File fh, void* buf, int count,
        MPI_Datatype datatype, MPI_Status* status);
int PMPI_File_iread_shared(MPI_File fh, void* buf, int count,
        MPI_Datatype datatype, MPI_Request* request);
int PMPI_File_iwrite_shared(MPI_File fh, void* buf, int count,
        MPI_Datatype datatype, MPI_Request* request);
int PMPI_File_read_ordered(MPI_File fh, void* buf, int count,
        MPI_Datatype datatype, MPI_Status* status);
int PMPI_File_write_ordered(MPI_File fh, void* buf, int count,
        MPI_Datatype datatype, MPI_Status* status);
int PMPI_File_seek_shared(MPI_File fh, MPI_Offset offset, int whence);
int PMPI_File_get_position_shared(MPI_File fh, MPI_Offset* offset);
int PMPI_File_read_at_all_begin(MPI_File fh, MPI_Offset offset, void* buf,
        int count, MPI_Datatype datatype);
int PMPI_File_read_at_all_end(MPI_File fh, void* buf, MPI_Status* status);
int PMPI_File_write_at_all_begin(MPI_File fh, MPI_Offset offset, void* buf,
        int count, MPI_Datatype datatype);
int PMPI_File_write_at_all_end(MPI_File fh, void* buf, MPI_Status* status);
int PMPI_File_read_all_begin(MPI_File fh, void* buf, int count,
        MPI_Datatype datatype);
int PMPI_File_read_all_end(MPI_File fh, void* buf, MPI_Status* status);
int PMPI_File_write_all_begin(MPI_File fh, void* buf, int count,
        MPI_Datatype datatype);
int PMPI_File_write_all_end(MPI_File fh, void* buf, MPI_Status* status);
int PMPI_File_read_ordered_begin(MPI_File fh, void* buf, int count,
        MPI_Datatype datatype);
int PMPI_File_read_ordered_end(MPI_File fh, void* buf, MPI_Status* status);
int PMPI_File_write_ordered_begin(MPI_File fh, void* buf, int count,
        MPI_Datatype datatype);
int PMPI_File_write_ordered_end(MPI_File fh, void* buf, MPI_Status* status);
int PMPI_File_get_type_extent(MPI_File fh, MPI_Datatype datatype,
        MPI_Aint* extent);
int PMPI_File_set_atomicity(MPI_File fh, int flag);
int PMPI_File_get_atomicity(MPI_File fh, int* flag);
int PMPI_File_sync(MPI_File fh);
int PMPI_Finalize();
int PMPI_Finalized(int* flag);
int PMPI_Free_mem(void* base);
int PMPI_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount, MPI_Datatype recvtype,
        int root, MPI_Comm comm);
int PMPI_Igather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount, MPI_Datatype recvtype,
        int root, MPI_Comm comm, MPI_Request* request);
int PMPI_Gatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int* recvcounts, int* displs,
        MPI_Datatype recvtype, int root, MPI_Comm comm);
int PMPI_Igatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int* recvcounts, int* displs,
        MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request* request);
int PMPI_Get_address(void* location, MPI_Aint* address);
int PMPI_Get_count(MPI_Status* status, MPI_Datatype datatype, int* count);
int PMPI_Get_elements(MPI_Status* status, MPI_Datatype datatype,
        int* count);
int PMPI_Get_elements_x(MPI_Status* status, MPI_Datatype datatype,
        MPI_Count* count);
int PMPI_Get(void* origin_addr, int origin_count,
        MPI_Datatype origin_datatype, int target_rank,
        MPI_Aint target_disp, int target_count,
        MPI_Datatype target_datatype, MPI_Win win);
int PMPI_Get_accumulate(void* origin_addr, int origin_count, MPI_Datatype origin_datatype,
        void* result_addr, int result_count, MPI_Datatype result_datatype,
        int target_rank, MPI_Aint target_disp, int target_count,
        MPI_Datatype target_datatype, MPI_Op op, MPI_Win win);
int PMPI_Get_library_version(char* version_, int* resultlen);
int PMPI_Get_processor_name(char* name, int* resultlen);
int PMPI_Get_version(int* version_, int* subversion);
int PMPI_Graph_create(MPI_Comm comm_old, int nnodes, int* index,
        int* edges, int reorder, MPI_Comm* comm_graph);
int PMPI_Graph_get(MPI_Comm comm, int maxindex, int maxedges,
        int* index, int* edges);
int PMPI_Graph_map(MPI_Comm comm, int nnodes, int* index, int* edges,
        int* newrank);
int PMPI_Graph_neighbors_count(MPI_Comm comm, int rank, int* nneighbors);
int PMPI_Graph_neighbors(MPI_Comm comm, int rank, int maxneighbors,
        int* neighbors);
int PMPI_Graphdims_get(MPI_Comm comm, int* nnodes, int* nedges);
int PMPI_Grequest_complete(MPI_Request request);
int PMPI_Grequest_start(MPI_Grequest_query_function* query_fn,
        MPI_Grequest_free_function* free_fn,
        MPI_Grequest_cancel_function* cancel_fn,
        void* extra_state, MPI_Request* request);
MPI_Fint PMPI_Group_c2f(MPI_Group group);
int PMPI_Group_compare(MPI_Group group1, MPI_Group group2, int* result);
int PMPI_Group_difference(MPI_Group group1, MPI_Group group2,
        MPI_Group* newgroup);
int PMPI_Group_excl(MPI_Group group, int n, int* ranks,
        MPI_Group* newgroup);
MPI_Group PMPI_Group_f2c(MPI_Fint group);
int PMPI_Group_free(MPI_Group* group);
int PMPI_Group_incl(MPI_Group group, int n, int* ranks,
        MPI_Group* newgroup);
int PMPI_Group_intersection(MPI_Group group1, MPI_Group group2,
        MPI_Group* newgroup);
int PMPI_Group_range_excl(MPI_Group group, int n, int[3]* ranges,
        MPI_Group* newgroup);
int PMPI_Group_range_incl(MPI_Group group, int n, int[3]* ranges,
        MPI_Group* newgroup);
int PMPI_Group_rank(MPI_Group group, int* rank);
int PMPI_Group_size(MPI_Group group, int* size);
int PMPI_Group_translate_ranks(MPI_Group group1, int n, int* ranks1,
        MPI_Group group2, int* ranks2);
int PMPI_Group_union(MPI_Group group1, MPI_Group group2,
        MPI_Group* newgroup);
int PMPI_Ibsend(void* buf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm, MPI_Request* request);
int PMPI_Improbe(int source, int tag, MPI_Comm comm,
        int* flag, MPI_Message* message,
        MPI_Status* status);
int PMPI_Imrecv(void* buf, int count, MPI_Datatype type,
        MPI_Message* message, MPI_Request* request);
MPI_Fint PMPI_Info_c2f(MPI_Info info);
int PMPI_Info_create(MPI_Info* info);
int PMPI_Info_delete(MPI_Info info, char* key);
int PMPI_Info_dup(MPI_Info info, MPI_Info* newinfo);
MPI_Info PMPI_Info_f2c(MPI_Fint info);
int PMPI_Info_free(MPI_Info* info);
int PMPI_Info_get(MPI_Info info, char* key, int valuelen,
        char* value, int* flag);
int PMPI_Info_get_nkeys(MPI_Info info, int* nkeys);
int PMPI_Info_get_nthkey(MPI_Info info, int n, char* key);
int PMPI_Info_get_valuelen(MPI_Info info, char* key, int* valuelen,
        int* flag);
int PMPI_Info_set(MPI_Info info, char* key, char* value);
int PMPI_Init(int* argc, char*** argv);
int PMPI_Initialized(int* flag);
int PMPI_Init_thread(int* argc, char*** argv, int required,
        int* provided);
int PMPI_Intercomm_create(MPI_Comm local_comm, int local_leader,
        MPI_Comm bridge_comm, int remote_leader,
        int tag, MPI_Comm* newintercomm);
int PMPI_Intercomm_merge(MPI_Comm intercomm, int high,
        MPI_Comm* newintercomm);
int PMPI_Iprobe(int source, int tag, MPI_Comm comm, int* flag,
        MPI_Status* status);
int PMPI_Irecv(void* buf, int count, MPI_Datatype datatype, int source,
        int tag, MPI_Comm comm, MPI_Request* request);
int PMPI_Irsend(void* buf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm, MPI_Request* request);
int PMPI_Isend(void* buf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm, MPI_Request* request);
int PMPI_Issend(void* buf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm, MPI_Request* request);
int PMPI_Is_thread_main(int* flag);
int PMPI_Lookup_name(char* service_name, MPI_Info info, char* port_name);
MPI_Fint PMPI_Message_c2f(MPI_Message message);
MPI_Message PMPI_Message_f2c(MPI_Fint message);
int PMPI_Mprobe(int source, int tag, MPI_Comm comm,
        MPI_Message* message,
        MPI_Status* status);
int PMPI_Mrecv(void* buf, int count, MPI_Datatype type,
        MPI_Message* message, MPI_Status* status);
int PMPI_Neighbor_allgather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount, MPI_Datatype recvtype,
        MPI_Comm comm);
int PMPI_Ineighbor_allgather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount, MPI_Datatype recvtype,
        MPI_Comm comm, MPI_Request* request);
int PMPI_Neighbor_allgatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int* recvcounts, int* displs,
        MPI_Datatype recvtype, MPI_Comm comm);
int PMPI_Ineighbor_allgatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int* recvcounts, int* displs,
        MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request);
int PMPI_Neighbor_alltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount, MPI_Datatype recvtype,
        MPI_Comm comm);
int PMPI_Ineighbor_alltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount, MPI_Datatype recvtype,
        MPI_Comm comm, MPI_Request* request);
int PMPI_Neighbor_alltoallv(void* sendbuf, int* sendcounts, int* sdispls,  MPI_Datatype sendtype,
        void* recvbuf, int* recvcounts, int* rdispls, MPI_Datatype recvtype,
        MPI_Comm comm);
int PMPI_Ineighbor_alltoallv(void* sendbuf, int* sendcounts, int* sdispls, MPI_Datatype sendtype,
        void* recvbuf, int* recvcounts, int* rdispls, MPI_Datatype recvtype,
        MPI_Comm comm, MPI_Request* request);
int PMPI_Neighbor_alltoallw(void* sendbuf, int* sendcounts, MPI_Aint* sdispls, MPI_Datatype* sendtypes,
        void* recvbuf, int* recvcounts, MPI_Aint* rdispls, MPI_Datatype* recvtypes,
        MPI_Comm comm);
int PMPI_Ineighbor_alltoallw(void* sendbuf, int* sendcounts, MPI_Aint* sdispls, MPI_Datatype* sendtypes,
        void* recvbuf, int* recvcounts, MPI_Aint* rdispls, MPI_Datatype* recvtypes,
        MPI_Comm comm, MPI_Request* request);
MPI_Fint PMPI_Op_c2f(MPI_Op op);
int PMPI_Op_commutative(MPI_Op op, int* commute);
int PMPI_Op_create(MPI_User_function* function_, int commute, MPI_Op* op);
int PMPI_Open_port(MPI_Info info, char* port_name);
MPI_Op PMPI_Op_f2c(MPI_Fint op);
int PMPI_Op_free(MPI_Op* op);
int PMPI_Pack_external(char* datarep, void* inbuf, int incount,
        MPI_Datatype datatype, void* outbuf,
        MPI_Aint outsize, MPI_Aint* position);
int PMPI_Pack_external_size(char* datarep, int incount,
        MPI_Datatype datatype, MPI_Aint* size);
int PMPI_Pack(void* inbuf, int incount, MPI_Datatype datatype,
        void* outbuf, int outsize, int* position, MPI_Comm comm);
int PMPI_Pack_size(int incount, MPI_Datatype datatype, MPI_Comm comm,
        int* size);
int PMPI_Pcontrol(int level, ...);
int PMPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status* status);
int PMPI_Publish_name(char* service_name, MPI_Info info,
        char* port_name);
int PMPI_Put(void* origin_addr, int origin_count, MPI_Datatype origin_datatype,
        int target_rank, MPI_Aint target_disp, int target_count,
        MPI_Datatype target_datatype, MPI_Win win);
int PMPI_Query_thread(int* provided);
int PMPI_Raccumulate(void* origin_addr, int origin_count, MPI_Datatype origin_datatype, 
        int target_rank, MPI_Aint target_disp, int target_count, 
        MPI_Datatype target_datatype, MPI_Op op, MPI_Win win, MPI_Request* request);
int PMPI_Recv_init(void* buf, int count, MPI_Datatype datatype, int source,
        int tag, MPI_Comm comm, MPI_Request* request);
int PMPI_Recv(void* buf, int count, MPI_Datatype datatype, int source,
        int tag, MPI_Comm comm, MPI_Status* status);
int PMPI_Reduce(void* sendbuf, void* recvbuf, int count,
        MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
int PMPI_Ireduce(void* sendbuf, void* recvbuf, int count,
        MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request* request);
int PMPI_Reduce_local(void* inbuf, void* inoutbuf, int count,
        MPI_Datatype datatype, MPI_Op);
int PMPI_Reduce_scatter(void* sendbuf, void* recvbuf, int* recvcounts,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int PMPI_Ireduce_scatter(void* sendbuf, void* recvbuf, int* recvcounts,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request* request);
int PMPI_Reduce_scatter_block(void* sendbuf, void* recvbuf, int recvcount,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int PMPI_Ireduce_scatter_block(void* sendbuf, void* recvbuf, int recvcount,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request* request);
int PMPI_Register_datarep(char* datarep,
        MPI_Datarep_conversion_function* read_conversion_fn,
        MPI_Datarep_conversion_function* write_conversion_fn,
        MPI_Datarep_extent_function* dtype_file_extent_fn,
        void* extra_state);
MPI_Fint PMPI_Request_c2f(MPI_Request request);
MPI_Request PMPI_Request_f2c(MPI_Fint request);
int PMPI_Request_free(MPI_Request* request);
int PMPI_Request_get_status(MPI_Request request, int* flag,
        MPI_Status* status);
int PMPI_Rget(void* origin_addr, int origin_count, MPI_Datatype origin_datatype, 
        int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype,
        MPI_Win win, MPI_Request* request);
int PMPI_Rget_accumulate(void* origin_addr, int origin_count, MPI_Datatype origin_datatype,
        void* result_addr, int result_count, MPI_Datatype result_datatype,
        int target_rank, MPI_Aint target_disp, int target_count, 
        MPI_Datatype target_datatype, MPI_Op op,
        MPI_Win win, MPI_Request* request);
int PMPI_Rput(void* origin_addr, int origin_count, MPI_Datatype origin_datatype,
        int target_rank, MPI_Aint target_disp, int target_cout, 
        MPI_Datatype target_datatype, MPI_Win win, MPI_Request* request);
int PMPI_Rsend(void* ibuf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm);
int PMPI_Rsend_init(void* buf, int count, MPI_Datatype datatype,
        int dest, int tag, MPI_Comm comm,
        MPI_Request* request);
int PMPI_Scan(void* sendbuf, void* recvbuf, int count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int PMPI_Iscan(void* sendbuf, void* recvbuf, int count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request* request);
int PMPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount, MPI_Datatype recvtype,
        int root, MPI_Comm comm);
int PMPI_Iscatter(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        void* recvbuf, int recvcount, MPI_Datatype recvtype,
        int root, MPI_Comm comm, MPI_Request* request);
int PMPI_Scatterv(void* sendbuf, int* sendcounts, int* displs,
        MPI_Datatype sendtype, void* recvbuf, int recvcount,
        MPI_Datatype recvtype, int root, MPI_Comm comm);
int PMPI_Iscatterv(void* sendbuf, int* sendcounts, int* displs,
        MPI_Datatype sendtype, void* recvbuf, int recvcount,
        MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request* request);
int PMPI_Send_init(void* buf, int count, MPI_Datatype datatype,
        int dest, int tag, MPI_Comm comm,
        MPI_Request* request);
int PMPI_Send(void* buf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm);
int PMPI_Sendrecv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
        int dest, int sendtag, void* recvbuf, int recvcount,
        MPI_Datatype recvtype, int source, int recvtag,
        MPI_Comm comm,  MPI_Status* status);
int PMPI_Sendrecv_replace(void*  buf, int count, MPI_Datatype datatype,
        int dest, int sendtag, int source, int recvtag,
        MPI_Comm comm, MPI_Status* status);
int PMPI_Ssend_init(void* buf, int count, MPI_Datatype datatype,
        int dest, int tag, MPI_Comm comm,
        MPI_Request* request);
int PMPI_Ssend(void* buf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm);
int PMPI_Start(MPI_Request* request);
int PMPI_Startall(int count, MPI_Request* array_of_requests);
int PMPI_Status_c2f(MPI_Status* c_status, MPI_Fint* f_status);
int PMPI_Status_f2c(MPI_Fint* f_status, MPI_Status* c_status);
int PMPI_Status_set_cancelled(MPI_Status* status, int flag);
int PMPI_Status_set_elements(MPI_Status* status, MPI_Datatype datatype,
        int count);
int PMPI_Status_set_elements_x(MPI_Status* status, MPI_Datatype datatype,
        MPI_Count count);
int PMPI_Testall(int count, MPI_Request* array_of_requests, int* flag,
        MPI_Status* array_of_statuses);
int PMPI_Testany(int count, MPI_Request* array_of_requests, int* index, int* flag, MPI_Status* status);
int PMPI_Test(MPI_Request* request, int* flag, MPI_Status* status);
int PMPI_Test_cancelled(MPI_Status* status, int* flag);
int PMPI_Testsome(int incount, MPI_Request* array_of_requests,
        int* outcount, int* array_of_indices,
        MPI_Status* array_of_statuses);
int PMPI_Topo_test(MPI_Comm comm, int* status);
MPI_Fint PMPI_Type_c2f(MPI_Datatype datatype);
int PMPI_Type_commit(MPI_Datatype* type);
int PMPI_Type_contiguous(int count, MPI_Datatype oldtype,
        MPI_Datatype* newtype);
int PMPI_Type_create_darray(int size, int rank, int ndims,
        int* gsize_array, int* distrib_array,
        int* darg_array, int* psize_array,
        int order, MPI_Datatype oldtype,
        MPI_Datatype* newtype);
int PMPI_Type_create_f90_complex(int p, int r, MPI_Datatype* newtype);
int PMPI_Type_create_f90_integer(int r, MPI_Datatype* newtype);
int PMPI_Type_create_f90_real(int p, int r, MPI_Datatype* newtype);
int PMPI_Type_create_hindexed(int count, int* array_of_blocklengths,
        MPI_Aint* array_of_displacements,
        MPI_Datatype oldtype,
        MPI_Datatype* newtype);
int PMPI_Type_create_hvector(int count, int blocklength, MPI_Aint stride,
        MPI_Datatype oldtype,
        MPI_Datatype* newtype);
int PMPI_Type_create_keyval(MPI_Type_copy_attr_function* type_copy_attr_fn,
        MPI_Type_delete_attr_function* type_delete_attr_fn,
        int* type_keyval, void* extra_state);
int PMPI_Type_create_hindexed_block(int count, int blocklength,
        MPI_Aint* array_of_displacements,
        MPI_Datatype oldtype,
        MPI_Datatype* newtype);
int PMPI_Type_create_indexed_block(int count, int blocklength,
        int* array_of_displacements,
        MPI_Datatype oldtype,
        MPI_Datatype* newtype);
int PMPI_Type_create_struct(int count, int* array_of_block_lengths,
        MPI_Aint* array_of_displacements,
        MPI_Datatype* array_of_types,
        MPI_Datatype* newtype);
int PMPI_Type_create_subarray(int ndims, int* size_array, int* subsize_array,
        int* start_array, int order,
        MPI_Datatype oldtype, MPI_Datatype* newtype);
int PMPI_Type_create_resized(MPI_Datatype oldtype, MPI_Aint lb,
        MPI_Aint extent, MPI_Datatype* newtype);
int PMPI_Type_delete_attr(MPI_Datatype type, int type_keyval);
int PMPI_Type_dup(MPI_Datatype type, MPI_Datatype* newtype);
int PMPI_Type_free(MPI_Datatype* type);
int PMPI_Type_free_keyval(int* type_keyval);
MPI_Datatype PMPI_Type_f2c(MPI_Fint datatype);
int PMPI_Type_get_attr(MPI_Datatype type, int type_keyval,
        void* attribute_val, int* flag);
int PMPI_Type_get_contents(MPI_Datatype mtype, int max_integers,
        int max_addresses, int max_datatypes,
        int* array_of_integers,
        MPI_Aint* array_of_addresses,
        MPI_Datatype* array_of_datatypes);
int PMPI_Type_get_envelope(MPI_Datatype type, int* num_integers,
        int* num_addresses, int* num_datatypes,
        int* combiner);
int PMPI_Type_get_extent(MPI_Datatype type, MPI_Aint* lb,
        MPI_Aint* extent);
int PMPI_Type_get_extent_x(MPI_Datatype type, MPI_Count* lb,
        MPI_Count* extent);
int PMPI_Type_get_name(MPI_Datatype type, char* type_name,
        int* resultlen);
int PMPI_Type_get_true_extent(MPI_Datatype datatype, MPI_Aint* true_lb,
        MPI_Aint* true_extent);
int PMPI_Type_get_true_extent_x(MPI_Datatype datatype, MPI_Count* true_lb,
        MPI_Count* true_extent);
int PMPI_Type_indexed(int count, int* array_of_blocklengths,
        int* array_of_displacements,
        MPI_Datatype oldtype, MPI_Datatype* newtype);
int PMPI_Type_match_size(int typeclass, int size, MPI_Datatype* type);
int PMPI_Type_set_attr(MPI_Datatype type, int type_keyval,
        void* attr_val);
int PMPI_Type_set_name(MPI_Datatype type, char* type_name);
int PMPI_Type_size(MPI_Datatype type, int* size);
int PMPI_Type_size_x(MPI_Datatype type, MPI_Count* size);
int PMPI_Type_vector(int count, int blocklength, int stride,
        MPI_Datatype oldtype, MPI_Datatype* newtype);
int PMPI_Unpack(void* inbuf, int insize, int* position,
        void* outbuf, int outcount, MPI_Datatype datatype,
        MPI_Comm comm);
int PMPI_Unpublish_name(char* service_name, MPI_Info info,
        char* port_name);
int PMPI_Unpack_external (char* datarep, void* inbuf, MPI_Aint insize,
        MPI_Aint* position, void* outbuf, int outcount,
        MPI_Datatype datatype);
int PMPI_Waitall(int count, MPI_Request* array_of_requests,
        MPI_Status* array_of_statuses);
int PMPI_Waitany(int count, MPI_Request* array_of_requests,
        int* index, MPI_Status* status);
int PMPI_Wait(MPI_Request* request, MPI_Status* status);
int PMPI_Waitsome(int incount, MPI_Request* array_of_requests,
        int* outcount, int* array_of_indices,
        MPI_Status* array_of_statuses);
int PMPI_Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info,
        MPI_Comm comm, void* baseptr, MPI_Win* win);
int PMPI_Win_allocate_shared(MPI_Aint size, int disp_unit, MPI_Info info,
        MPI_Comm comm, void* baseptr, MPI_Win* win);
int PMPI_Win_attach(MPI_Win win, void* base, MPI_Aint size);
MPI_Fint PMPI_Win_c2f(MPI_Win win);
int PMPI_Win_call_errhandler(MPI_Win win, int errorcode);
int PMPI_Win_complete(MPI_Win win);
int PMPI_Win_create(void* base, MPI_Aint size, int disp_unit,
        MPI_Info info, MPI_Comm comm, MPI_Win* win);
int PMPI_Win_create_dynamic(MPI_Info info, MPI_Comm comm, MPI_Win* win);
int PMPI_Win_create_errhandler(MPI_Win_errhandler_function* function_,
        MPI_Errhandler* errhandler);
int PMPI_Win_create_keyval(MPI_Win_copy_attr_function* win_copy_attr_fn,
        MPI_Win_delete_attr_function* win_delete_attr_fn,
        int* win_keyval, void* extra_state);
int PMPI_Win_delete_attr(MPI_Win win, int win_keyval);
int PMPI_Win_detach(MPI_Win win, void* base);
MPI_Win PMPI_Win_f2c(MPI_Fint win);
int PMPI_Win_fence(int assert_, MPI_Win win);
int PMPI_Win_flush(int rank, MPI_Win win);
int PMPI_Win_flush_all(MPI_Win win);
int PMPI_Win_flush_local(int rank, MPI_Win win);
int PMPI_Win_flush_local_all(MPI_Win win);
int PMPI_Win_free(MPI_Win* win);
int PMPI_Win_free_keyval(int* win_keyval);
int PMPI_Win_get_attr(MPI_Win win, int win_keyval,
        void* attribute_val, int* flag);
int PMPI_Win_get_errhandler(MPI_Win win, MPI_Errhandler* errhandler);
int PMPI_Win_get_group(MPI_Win win, MPI_Group* group);
int PMPI_Win_get_info(MPI_Win win, MPI_Info* info_used);
int PMPI_Win_get_name(MPI_Win win, char* win_name, int* resultlen);
int PMPI_Win_lock(int lock_type, int rank, int assert_, MPI_Win win);
int PMPI_Win_lock_all(int assert_, MPI_Win win);
int PMPI_Win_post(MPI_Group group, int assert_, MPI_Win win);
int PMPI_Win_set_attr(MPI_Win win, int win_keyval, void* attribute_val);
int PMPI_Win_set_errhandler(MPI_Win win, MPI_Errhandler errhandler);
int PMPI_Win_set_info(MPI_Win win, MPI_Info info);
int PMPI_Win_set_name(MPI_Win win, char* win_name);
int PMPI_Win_shared_query(MPI_Win win, int rank, MPI_Aint* size, int* disp_unit, void* baseptr);
int PMPI_Win_start(MPI_Group group, int assert_, MPI_Win win);
int PMPI_Win_sync(MPI_Win win);
int PMPI_Win_test(MPI_Win win, int* flag);
int PMPI_Win_unlock(int rank, MPI_Win win);
int PMPI_Win_unlock_all(MPI_Win win);
int PMPI_Win_wait(MPI_Win win);
double PMPI_Wtick();
double PMPI_Wtime();
int PMPI_T_init_thread (int required, int* provided);
int PMPI_T_finalize ();
int PMPI_T_cvar_get_num (int* num_cvar);
int PMPI_T_cvar_get_info (int cvar_index, char* name, int* name_len,
        int* verbosity, MPI_Datatype* datatype,
        MPI_T_enum* enumtype, char* desc,
        int* desc_len, int* bind, int* scope_);
int PMPI_T_cvar_get_index (char* name, int* cvar_index);
int PMPI_T_cvar_handle_alloc (int cvar_index, void* obj_handle,
        MPI_T_cvar_handle* handle, int* count);
int PMPI_T_cvar_handle_free (MPI_T_cvar_handle* handle);
int PMPI_T_cvar_read (MPI_T_cvar_handle handle, void* buf);
int PMPI_T_cvar_write (MPI_T_cvar_handle handle, void* buf);
int PMPI_T_category_get_num(int* num_cat);
int PMPI_T_category_get_info(int cat_index, char* name, int* name_len,
        char* desc, int* desc_len, int* num_cvars,
        int* num_pvars, int* num_categories);
int PMPI_T_category_get_index (char* name, int* category_index);
int PMPI_T_category_get_cvars(int cat_index, int len, int* indices);
int PMPI_T_category_get_pvars(int cat_index, int len, int* indices);
int PMPI_T_category_get_categories(int cat_index, int len, int* indices);
int PMPI_T_category_changed(int* stamp);
    
int PMPI_T_pvar_get_num(int* num_pvar);
int PMPI_T_pvar_get_info(int pvar_index, char* name, int* name_len,
        int* verbosity, int* var_class, MPI_Datatype* datatype,
        MPI_T_enum* enumtype, char* desc, int* desc_len, int* bind,
        int* readonly, int* continuous, int* atomic);
int PMPI_T_pvar_get_index (char* name, int var_class, int* pvar_index);
int PMPI_T_pvar_session_create(MPI_T_pvar_session* session);
int PMPI_T_pvar_session_free(MPI_T_pvar_session* session);
int PMPI_T_pvar_handle_alloc(MPI_T_pvar_session session, int pvar_index,
        void* obj_handle, MPI_T_pvar_handle* handle, int* count);
int PMPI_T_pvar_handle_free(MPI_T_pvar_session session, MPI_T_pvar_handle* handle);
int PMPI_T_pvar_start(MPI_T_pvar_session session, MPI_T_pvar_handle handle);
int PMPI_T_pvar_stop(MPI_T_pvar_session session, MPI_T_pvar_handle handle);
int PMPI_T_pvar_read(MPI_T_pvar_session session, MPI_T_pvar_handle handle,
        void* buf);
int PMPI_T_pvar_write(MPI_T_pvar_session session, MPI_T_pvar_handle handle,
        void* buf);
int PMPI_T_pvar_reset(MPI_T_pvar_session session, MPI_T_pvar_handle handle);
int PMPI_T_pvar_readreset(MPI_T_pvar_session session, MPI_T_pvar_handle handle,
        void* buf);
int PMPI_T_enum_get_info(MPI_T_enum enumtype, int* num, char* name, int* name_len);
int PMPI_T_enum_get_item(MPI_T_enum enumtype, int index, int* value, char* name,
        int* name_len);

/*
 *  Tool MPI API
 */
int MPI_T_init_thread (int required, int* provided);
int MPI_T_finalize ();
int MPI_T_cvar_get_num (int* num_cvar);
int MPI_T_cvar_get_info (int cvar_index, char* name, int* name_len,
        int* verbosity, MPI_Datatype* datatype,
        MPI_T_enum* enumtype, char* desc,
        int* desc_len, int* bind, int* scope_);
int MPI_T_cvar_get_index (char* name, int* cvar_index);
int MPI_T_cvar_handle_alloc (int cvar_index, void* obj_handle,
        MPI_T_cvar_handle* handle, int* count);
int MPI_T_cvar_handle_free (MPI_T_cvar_handle* handle);
int MPI_T_cvar_read (MPI_T_cvar_handle handle, void* buf);
int MPI_T_cvar_write (MPI_T_cvar_handle handle, void* buf);
int MPI_T_category_get_num(int* num_cat);
int MPI_T_category_get_info(int cat_index, char* name, int* name_len,
        char* desc, int* desc_len, int* num_cvars,
        int* num_pvars, int* num_categories);
int MPI_T_category_get_index (char* name, int* category_index);
int MPI_T_category_get_cvars(int cat_index, int len, int* indices);
int MPI_T_category_get_pvars(int cat_index, int len, int* indices);
int MPI_T_category_get_categories(int cat_index, int len, int* indices);
int MPI_T_category_changed(int* stamp);

int MPI_T_pvar_get_num(int* num_pvar);
int MPI_T_pvar_get_info(int pvar_index, char* name, int* name_len,
       int* verbosity, int* var_class, MPI_Datatype* datatype,
       MPI_T_enum* enumtype, char* desc, int* desc_len, int* bind,
       int* readonly, int* continuous, int* atomic);
int MPI_T_pvar_get_index (char* name, int var_class, int* pvar_index);
int MPI_T_pvar_session_create(MPI_T_pvar_session* session);
int MPI_T_pvar_session_free(MPI_T_pvar_session* session);
int MPI_T_pvar_handle_alloc(MPI_T_pvar_session session, int pvar_index,
        void* obj_handle, MPI_T_pvar_handle* handle, int* count);
int MPI_T_pvar_handle_free(MPI_T_pvar_session session, MPI_T_pvar_handle* handle);
int MPI_T_pvar_start(MPI_T_pvar_session session, MPI_T_pvar_handle handle);
int MPI_T_pvar_stop(MPI_T_pvar_session session, MPI_T_pvar_handle handle);
int MPI_T_pvar_read(MPI_T_pvar_session session, MPI_T_pvar_handle handle,
        void* buf);
int MPI_T_pvar_write(MPI_T_pvar_session session, MPI_T_pvar_handle handle,
        void* buf);
int MPI_T_pvar_reset(MPI_T_pvar_session session, MPI_T_pvar_handle handle);
int MPI_T_pvar_readreset(MPI_T_pvar_session session, MPI_T_pvar_handle handle,
        void* buf);
int MPI_T_enum_get_info(MPI_T_enum enumtype, int* num, char* name, int* name_len);
int MPI_T_enum_get_item(MPI_T_enum enumtype, int index, int* value, char* name,
        int* name_len);
