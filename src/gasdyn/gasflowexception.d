// gasflowexception.d

module gasflowexception;


class GasFlowException : Exception {
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}
