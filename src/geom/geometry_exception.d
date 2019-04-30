/**
 * Authors: Rowan G. and Peter J.
 * Date: 2019-04-30
 *
 */

module geom.geometry_exception;

class GeometryException : Exception {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}
