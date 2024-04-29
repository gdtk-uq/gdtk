/**
 * Module to hold Eilmer-specific exceptions.
 *
 * Exceptions are:
 *  + LmrException : a non-specific exception
 *  + NewtonKrylovExcepetion : when things go awry in the N-K solver
 *  + TimeMarchingException : for signalling problems in the time-marcher
 *  + LmrPostProcessingException : for post-processing badness
 *  + LmrPreProcessingException : for pre-processing badness
 *  + LmrInitialisationException : for bad events during the simulation (D-lang) initialisation
 *
 */

module lmrexceptions;

/**
 * Base Eilmer exception.
 *
 * This exception serves as a base class for all Eilmer exceptions.
 *
 * Author: RJG
 * Date: 2023-05-07
 */
class LmrException : Exception {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}

/**
 * Class to signal N-K specific exceptions.
 *
 * Author: RJG
 * Date: 2023-05-07
 */
class NewtonKrylovException : LmrException {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}

/**
 * Class to signal timemarching specific exceptions.
 *
 * Author: RJG
 * Date: 2024-03-09
 */
class TimeMarchingException : LmrException {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}

/**
 * Class to signal an exception at post-processing stage.
 *
 * Authors: RJG
 * Date: 2024-03-04
 */
class LmrPostProcessingException : LmrException {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}

/**
 * Class to signal an exception at pre-processing stage.
 *
 * Authors: RJG
 * Date: 2024-03-09
 */
class LmrPreProcessingException : LmrException {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}

/**
 * Class to signal an exception at run-time initialisation.
 *
 * Authors: RJG
 * Date: 2024-04-27
 */
class LmrInitialisationException : LmrException {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}
