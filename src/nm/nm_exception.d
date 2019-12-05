/**
 * An exception to signal problems specifically
 * in the numerical methods functions.
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2019-12-05
 */

module nm.nm_exception;

class NumericalMethodException : Exception {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}
