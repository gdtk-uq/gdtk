// gasflowexception.d

module gasdyn.gasflowexception;


class GasFlowException : Exception {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}
