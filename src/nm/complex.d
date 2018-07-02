// Written in the D programming language.

/** This module contains the $(LREF Complex) type, which is used to represent
    _complex numbers, along with related mathematical operations and functions.

    $(LREF Complex) will eventually
    $(DDLINK deprecate, Deprecated Features, replace)
    the built-in types `cfloat`, `cdouble`, `creal`, `ifloat`,
    `idouble`, and `ireal`.

    Authors:    Lars Tandle Kyllingstad, Don Clugston
    Copyright:  Copyright (c) 2010, Lars T. Kyllingstad.
    License:    $(HTTP boost.org/LICENSE_1_0.txt, Boost License 1.0)
    Source:     $(PHOBOSSRC std/_complex.d)

    Adapted for use in complexifying the Eilmer code
    Kyle Damm, 2018
*/
module nm.complex;

version(complex_numbers) {
    // PJ's 2018-05-28 hack to avoid compiler problems for most users.
    //
    // We include the body of this file, only if we really are building the
    // complex_numbers verion of the code because we wish to avoid
    // the internal compiler error that seems to happen for DMD v2.76+ compilers
    // when optimization is requested.
    // Look toward the end of this file for the brief section for the build
    // assuming double_numbers.  There we only need the simple structure definition.

/** Helper function that returns a _complex number with the specified
    real and imaginary parts.

    Params:
        R = (template parameter) type of real part of complex number
        I = (template parameter) type of imaginary part of complex number

        re = real part of complex number to be constructed
        im = (optional) imaginary part of complex number, 0 if omitted.

    Returns:
        `Complex` instance with real and imaginary parts set
        to the values provided as input.  If neither `re` nor
        `im` are floating-point numbers, the return type will
        be `Complex!double`.  Otherwise, the return type is
        deduced using $(D std.traits.CommonType!(R, I)).
*/

import std.traits;
import std.math;
import std.conv;

auto complex(R)(R re)  @safe pure nothrow @nogc
if (is(R : double))
{
    static if (isFloatingPoint!R)
        return Complex!R(re, 0);
    else
        return Complex!double(re, 0);
}

/// ditto
auto complex(R, I)(R re, I im)  @safe pure nothrow @nogc
if (is(R : double) && is(I : double))
{
    static if (isFloatingPoint!R || isFloatingPoint!I)
        return Complex!(CommonType!(R, I))(re, im);
    else
        return Complex!double(re, im);
}

///
@safe pure nothrow unittest
{
    auto a = complex(1.0);
    static assert(is(typeof(a) == Complex!double));
    assert(a.re == 1.0);
    assert(a.im == 0.0);

    auto b = complex(2.0L);
    static assert(is(typeof(b) == Complex!real));
    assert(b.re == 2.0L);
    assert(b.im == 0.0L);

    auto c = complex(1.0, 2.0);
    static assert(is(typeof(c) == Complex!double));
    assert(c.re == 1.0);
    assert(c.im == 2.0);

    auto d = complex(3.0, 4.0L);
    static assert(is(typeof(d) == Complex!real));
    assert(d.re == 3.0);
    assert(d.im == 4.0L);

    auto e = complex(1);
    static assert(is(typeof(e) == Complex!double));
    assert(e.re == 1);
    assert(e.im == 0);

    auto f = complex(1L, 2);
    static assert(is(typeof(f) == Complex!double));
    assert(f.re == 1L);
    assert(f.im == 2);

    auto g = complex(3, 4.0L);
    static assert(is(typeof(g) == Complex!real));
    assert(g.re == 3);
    assert(g.im == 4.0L);
}

/** A complex number parametrised by a type `T`, which must be either
    `float`, `double` or `real`.
*/
struct Complex(T)
if (isFloatingPoint!T)
{
    import std.format : FormatSpec;
    import std.range.primitives : isOutputRange;

    /** The real part of the number. */
    T re;

    /** The imaginary part of the number. */
    T im;

    /** Converts the complex number to a string representation.

    The second form of this function is usually not called directly;
    instead, it is used via $(REF format, std,string), as shown in the examples
    below.  Supported format characters are 'e', 'f', 'g', 'a', and 's'.

    See the $(MREF std, format) and $(REF format, std,string)
    documentation for more information.
    */
    string toString() const @safe /* TODO: pure nothrow */
    {
        import std.exception : assumeUnique;
        char[] buf;
        buf.reserve(100);
        auto fmt = FormatSpec!char("%s");
        toString((const(char)[] s) { buf ~= s; }, fmt);
        static trustedAssumeUnique(T)(T t) @trusted { return assumeUnique(t); }
        return trustedAssumeUnique(buf);
    }

    static if (is(T == double))
    ///
    @safe unittest
    {
        auto c = complex(1.2, 3.4);

        // Vanilla toString formatting:
        assert(c.toString() == "1.2+3.4i");

        // Formatting with std.string.format specs: the precision and width
        // specifiers apply to both the real and imaginary parts of the
        // complex number.
        import std.format : format;
        assert(format("%.2f", c)  == "1.20+3.40i");
        assert(format("%4.1f", c) == " 1.2+ 3.4i");
    }

    /// ditto
    void toString(Writer, Char)(scope Writer w, const ref FormatSpec!Char formatSpec) const
        if (isOutputRange!(Writer, const(Char)[]))
    {
        import std.format : formatValue;
        import std.math : signbit;
        import std.range.primitives : put;
        formatValue(w, re, formatSpec);
        if (signbit(im) == 0)
           put(w, "+");
        formatValue(w, im, formatSpec);
        put(w, "i");
    }

    // Construct from a string, assuming that we only ever want to set real part.
    // Note that to!T is not @nogc, so we put it before the other constructors.
    // PJ 2018-06-02
    this(S : string)(S s)
    {
        re = to!T(s);
        im = 0;
    }
    // Assignment operator
    // this = string
    ref Complex opAssign(S : string)(S s)
    {
        re = to!T(s);
        im = 0;
        return this;
    }

@safe pure nothrow @nogc:

    /** Construct a complex number with the specified real and
    imaginary parts. In the case where a single argument is passed
    that is not complex, the imaginary part of the result will be
    zero.
    */
    this(R : T)(Complex!R z)
    {
        re = z.re;
        im = z.im;
    }

    /// ditto
    this(Rx : T, Ry : T)(Rx x, Ry y)
    {
        re = x;
        im = y;
    }

    /// ditto
    this(R : T)(R r)
    {
        re = r;
        im = 0;
    }
    
    // ASSIGNMENT OPERATORS

    // this = complex
    ref Complex opAssign(R : T)(Complex!R z)
    {
        re = z.re;
        im = z.im;
        return this;
    }

    // this = numeric
    ref Complex opAssign(R : T)(R r)
    {
        re = r;
        im = 0;
        return this;
    }

    // COMPARISON OPERATORS

    // this == complex
    bool opEquals(R : T)(Complex!R z) const
    {
        return re == z.re;// && im == z.im;
    }

    // this == numeric
    bool opEquals(R : T)(R r) const
    {
        return re == r; // && im == 0;
    }

    // OpCmp structure from sources:
    // https://forum.dlang.org/post/p1qa9n$2hje$1@digitalmars.com
    // http://www.angelcode.com/angelscript/sdk/docs/manual/doc_script_class_ops.html

    // this = complex (KD, 2018)
    int opCmp(Complex!double z) const
    {
        auto diff = re - z.re;
        if (  std.math.fabs(diff) < 1.0e-50 )
            return 0;
        else if ( diff < 0 )
            return -1;
        else
            return 1;
    }

    // this = numeric (KD, 2018)
    int opCmp(double z) const
    {
        auto diff = re - z.re;
        if (  std.math.fabs(diff) < 1.0e-50 )
            return 0;
        else if ( diff < 0 )
            return -1;
        else
            return 1;
    }

    // UNARY OPERATORS

    // +complex
    Complex opUnary(string op)() const
        if (op == "+")
    {
        return this;
    }

    // -complex
    Complex opUnary(string op)() const
        if (op == "-")
    {
        return Complex(-re, -im);
    }

    // BINARY OPERATORS

    // complex op complex
    Complex!(CommonType!(T,R)) opBinary(string op, R)(Complex!R z) const
    {
        alias C = typeof(return);
        auto w = C(this.re, this.im);
        return w.opOpAssign!(op)(z);
    }

    // complex op numeric
    Complex!(CommonType!(T,R)) opBinary(string op, R)(R r) const
        if (isNumeric!R)
    {
        alias C = typeof(return);
        auto w = C(this.re, this.im);
        return w.opOpAssign!(op)(r);
    }

    // numeric + complex,  numeric * complex
    Complex!(CommonType!(T, R)) opBinaryRight(string op, R)(R r) const
        if ((op == "+" || op == "*") && (isNumeric!R))
    {
        return opBinary!(op)(r);
    }

    // numeric - complex
    Complex!(CommonType!(T, R)) opBinaryRight(string op, R)(R r) const
        if (op == "-" && isNumeric!R)
    {
        return Complex(r - re, -im);
    }

    // numeric / complex
    Complex!(CommonType!(T, R)) opBinaryRight(string op, R)(R r) const
        if (op == "/" && isNumeric!R)
    {
        import std.math : fabs;
        typeof(return) w = void;
        if (fabs(re) < fabs(im))
        {
            immutable ratio = re/im;
            immutable rdivd = r/(re*ratio + im);

            w.re = rdivd*ratio;
            w.im = -rdivd;
        }
        else
        {
            immutable ratio = im/re;
            immutable rdivd = r/(re + im*ratio);

            w.re = rdivd;
            w.im = -rdivd*ratio;
        }

        return w;
    }

    // numeric ^^ complex
    Complex!(CommonType!(T, R)) opBinaryRight(string op, R)(R lhs) const
        if (op == "^^" && isNumeric!R)
    {
        import std.math : cos, exp, log, sin, PI;
        Unqual!(CommonType!(T, R)) ab = void, ar = void;

        if (lhs >= 0)
        {
            // r = lhs
            // theta = 0
            ab = lhs ^^ this.re;
            ar = log(lhs) * this.im;
        }
        else
        {
            // r = -lhs
            // theta = PI
            ab = (-lhs) ^^ this.re * exp(-PI * this.im);
            ar = PI * this.re + log(-lhs) * this.im;
        }

        return typeof(return)(ab * cos(ar), ab * sin(ar));
    }

    // OP-ASSIGN OPERATORS

    // complex += complex,  complex -= complex
    ref Complex opOpAssign(string op, C)(C z)
        if ((op == "+" || op == "-") && is(C R == Complex!R))
    {
        mixin ("re "~op~"= z.re;");
        mixin ("im "~op~"= z.im;");
        return this;
    }

    // complex *= complex
    ref Complex opOpAssign(string op, C)(C z)
        if (op == "*" && is(C R == Complex!R))
    {
        auto temp = re*z.re - im*z.im;
        im = im*z.re + re*z.im;
        re = temp;
        return this;
    }

    // complex /= complex
    ref Complex opOpAssign(string op, C)(C z)
        if (op == "/" && is(C R == Complex!R))
    {
        import std.math : fabs;
        if (fabs(z.re) < fabs(z.im))
        {
            immutable ratio = z.re/z.im;
            immutable denom = z.re*ratio + z.im;

            immutable temp = (re*ratio + im)/denom;
            im = (im*ratio - re)/denom;
            re = temp;
        }
        else
        {
            immutable ratio = z.im/z.re;
            immutable denom = z.re + z.im*ratio;

            immutable temp = (re + im*ratio)/denom;
            im = (im - re*ratio)/denom;
            re = temp;
        }
        return this;
    }

    // complex ^^= complex
    // same as pow(Complex!double z, Complex!double w)
    ref Complex opOpAssign(string op, C)(C w)
        if (op == "^^" && is(C R == Complex!R))
    {
        double a = this.re; double b = this.im;
        double c = w.re; double d = w.im;
        double r = std.math.sqrt(a*a + b*b);
        double theta = arg(this); Complex!double i = complex(0.0, 1.0);
        double logr = std.math.log(to!double(r));
        re = exp(logr*w+i*theta*w).re; 
        im = exp(logr*w+i*theta*w).im;
        return this;
    }

    // removed (KD, 2018)
    /*
      ref Complex opOpAssign(string op, C)(C z)
        if (op == "^^" && is(C R == Complex!R))
    {
        import std.math : exp, log, cos, sin;
        immutable r = abs(this);
        immutable t = arg(this);
        immutable ab = r^^z.re * exp(-t*z.im);
        immutable ar = t*z.re + log(r)*z.im;

        re = ab*cos(ar);
        im = ab*sin(ar);
        return this;
    }
    */
    
    // complex += numeric,  complex -= numeric
    ref Complex opOpAssign(string op, U : T)(U a)
        if (op == "+" || op == "-")
    {
        mixin ("re "~op~"= a;");
        return this;
    }

    // complex *= numeric,  complex /= numeric
    ref Complex opOpAssign(string op, U : T)(U a)
        if (op == "*" || op == "/")
    {
        mixin ("re "~op~"= a;");
        mixin ("im "~op~"= a;");
        return this;
    }

    // complex ^^= real (KD, 2018)
    ref Complex opOpAssign(string op, R)(R r)
        if (op == "^^" && isFloatingPoint!R)
    {
        double a = this.re; double b = this.im;
        double p = std.math.sqrt(a*a + b*b);
        double logp = std.math.log(to!double(p));
        double theta = arg(this); Complex!double i = complex(0.0, 1.0);
        re = exp(logp*r+i*theta*r).re;
        im = exp(logp*r+i*theta*r).im;
        return this;
    }

    // removed (KD, 2018)
    /*
    ref Complex opOpAssign(string op, R)(R r)
        if (op == "^^" && isFloatingPoint!R)
    {
        import std.math : cos, sin;
        immutable ab = abs(this)^^r;
        immutable ar = arg(this)*r;
        re = ab*cos(ar);
        im = ab*sin(ar);
        return this;
    }
    */

    // complex ^^= int (KD, 2018)
    // same as pow(Complex!double z, int w)
    ref Complex opOpAssign(string op, U)(U l)
        if (op == "^^" && isIntegral!U)
    {
        Complex!double p  = complex(this.re, this.im);
        foreach ( i; 0..l-1) this *= p;
        return this;
    }

    // removed (KD, 2018)
    /*
    ref Complex opOpAssign(string op, U)(U i)
        if (op == "^^" && isIntegral!U)
    {
        switch (i)
        {
        case 0:
            re = 1.0;
            im = 0.0;
            break;
        case 1:
            // identity; do nothing
            break;
        case 2:
            this *= this;
            break;
        case 3:
            auto z = this;
            this *= z;
            this *= z;
            break;
        default:
            this ^^= cast(real) i;
        }
        return this;
    }
    */
}

@safe pure nothrow unittest
{
    import std.complex;
    import std.math;

    enum EPS = double.epsilon;
    auto c1 = complex(1.0, 1.0);

    // Check unary operations.
    auto c2 = Complex!double(0.5, 2.0);

    assert(c2 == +c2);

    assert((-c2).re == -(c2.re));
    assert((-c2).im == -(c2.im));
    assert(c2 == -(-c2));

    // Check complex-complex operations.
    auto cpc = c1 + c2;
    assert(cpc.re == c1.re + c2.re);
    assert(cpc.im == c1.im + c2.im);

    auto cmc = c1 - c2;
    assert(cmc.re == c1.re - c2.re);
    assert(cmc.im == c1.im - c2.im);

    auto ctc = c1 * c2;
    assert(approxEqual(abs(ctc), abs(c1)*abs(c2), EPS));
    assert(approxEqual(arg(ctc), arg(c1)+arg(c2), EPS));

    auto cdc = c1 / c2;
    assert(approxEqual(abs(cdc), abs(c1)/abs(c2), EPS));
    assert(approxEqual(arg(cdc), arg(c1)-arg(c2), EPS));

    auto cec = c1^^c2;
    assert(approxEqual(cec.re, 0.11524131979943839881, EPS));
    assert(approxEqual(cec.im, 0.21870790452746026696, EPS));

    // Check complex-real operations.
    double a = 123.456;

    auto cpr = c1 + a;
    assert(cpr.re == c1.re + a);
    assert(cpr.im == c1.im);

    auto cmr = c1 - a;
    assert(cmr.re == c1.re - a);
    assert(cmr.im == c1.im);

    auto ctr = c1 * a;
    assert(ctr.re == c1.re*a);
    assert(ctr.im == c1.im*a);

    auto cdr = c1 / a;
    assert(approxEqual(abs(cdr), abs(c1)/a, EPS));
    assert(approxEqual(arg(cdr), arg(c1), EPS));

    auto cer = c1^^3.0;
    assert(approxEqual(abs(cer), abs(c1)^^3, EPS));
    assert(approxEqual(arg(cer), arg(c1)*3, EPS));

    auto rpc = a + c1;
    assert(rpc == cpr);

    auto rmc = a - c1;
    assert(rmc.re == a-c1.re);
    assert(rmc.im == -c1.im);

    auto rtc = a * c1;
    assert(rtc == ctr);

    auto rdc = a / c1;
    assert(approxEqual(abs(rdc), a/abs(c1), EPS));
    assert(approxEqual(arg(rdc), -arg(c1), EPS));

    rdc = a / c2;
    assert(approxEqual(abs(rdc), a/abs(c2), EPS));
    assert(approxEqual(arg(rdc), -arg(c2), EPS));

    auto rec1a = 1.0 ^^ c1;
    assert(rec1a.re == 1.0);
    assert(rec1a.im == 0.0);

    auto rec2a = 1.0 ^^ c2;
    assert(rec2a.re == 1.0);
    assert(rec2a.im == 0.0);

    auto rec1b = (-1.0) ^^ c1;
    assert(approxEqual(abs(rec1b), std.math.exp(-PI * c1.im), EPS));
    auto arg1b = arg(rec1b);
    /* The argument _should_ be PI, but floating-point rounding error
     * means that in fact the imaginary part is very slightly negative.
     */
    assert(approxEqual(arg1b, PI, EPS) || approxEqual(arg1b, -PI, EPS));

    auto rec2b = (-1.0) ^^ c2;
    assert(approxEqual(abs(rec2b), std.math.exp(-2 * PI), EPS));
    assert(approxEqual(arg(rec2b), PI_2, EPS));

    auto rec3a = 0.79 ^^ complex(6.8, 5.7);
    auto rec3b = complex(0.79, 0.0) ^^ complex(6.8, 5.7);
    assert(approxEqual(rec3a.re, rec3b.re, EPS));
    assert(approxEqual(rec3a.im, rec3b.im, EPS));

    auto rec4a = (-0.79) ^^ complex(6.8, 5.7);
    auto rec4b = complex(-0.79, 0.0) ^^ complex(6.8, 5.7);
    assert(approxEqual(rec4a.re, rec4b.re, EPS));
    assert(approxEqual(rec4a.im, rec4b.im, EPS));

    auto rer = a ^^ complex(2.0, 0.0);
    auto rcheck = a ^^ 2.0;
    static assert(is(typeof(rcheck) == double));
    assert(feqrel(rer.re, rcheck) == double.mant_dig);
    assert(isIdentical(rer.re, rcheck));
    assert(rer.im == 0.0);

    auto rer2 = (-a) ^^ complex(2.0, 0.0);
    rcheck = (-a) ^^ 2.0;
    assert(feqrel(rer2.re, rcheck) == double.mant_dig);
    assert(isIdentical(rer2.re, rcheck));
    assert(approxEqual(rer2.im, 0.0, EPS));

    auto rer3 = (-a) ^^ complex(-2.0, 0.0);
    rcheck = (-a) ^^ (-2.0);
    assert(feqrel(rer3.re, rcheck) == double.mant_dig);
    assert(isIdentical(rer3.re, rcheck));
    assert(approxEqual(rer3.im, 0.0, EPS));

    auto rer4 = a ^^ complex(-2.0, 0.0);
    rcheck = a ^^ (-2.0);
    assert(feqrel(rer4.re, rcheck) == double.mant_dig);
    assert(isIdentical(rer4.re, rcheck));
    assert(rer4.im == 0.0);

    // Check Complex-int operations.
    foreach (i; 0 .. 6)
    {
        auto cei = c1^^i;
        assert(approxEqual(abs(cei), abs(c1)^^i, EPS));
        // Use cos() here to deal with arguments that go outside
        // the (-pi,pi] interval (only an issue for i>3).
        assert(approxEqual(std.math.cos(arg(cei)), std.math.cos(arg(c1)*i), EPS));
    }

    // Check operations between different complex types.
    auto cf = Complex!float(1.0, 1.0);
    auto cr = Complex!real(1.0, 1.0);
    auto c1pcf = c1 + cf;
    auto c1pcr = c1 + cr;
    static assert(is(typeof(c1pcf) == Complex!double));
    static assert(is(typeof(c1pcr) == Complex!real));
    assert(c1pcf.re == c1pcr.re);
    assert(c1pcf.im == c1pcr.im);

    auto c1c = c1;
    auto c2c = c2;

    c1c /= c1;
    assert(approxEqual(c1c.re, 1.0, EPS));
    assert(approxEqual(c1c.im, 0.0, EPS));

    c1c = c1;
    c1c /= c2;
    assert(approxEqual(c1c.re, 0.588235, EPS));
    assert(approxEqual(c1c.im, -0.352941, EPS));

    c2c /= c1;
    assert(approxEqual(c2c.re, 1.25, EPS));
    assert(approxEqual(c2c.im, 0.75, EPS));

    c2c = c2;
    c2c /= c2;
    assert(approxEqual(c2c.re, 1.0, EPS));
    assert(approxEqual(c2c.im, 0.0, EPS));
}

@safe pure nothrow unittest
{
    // Initialization
    Complex!double a = 1;
    assert(a.re == 1 && a.im == 0);
    Complex!double b = 1.0;
    assert(b.re == 1.0 && b.im == 0);
    Complex!double c = Complex!real(1.0, 2);
    assert(c.re == 1.0 && c.im == 2);
}

@safe pure nothrow unittest
{
    // Assignments and comparisons
    Complex!double z;

    z = 1;
    assert(z == 1);
    assert(z.re == 1.0  &&  z.im == 0.0);

    z = 2.0;
    assert(z == 2.0);
    assert(z.re == 2.0  &&  z.im == 0.0);

    z = 1.0L;
    assert(z == 1.0L);
    assert(z.re == 1.0  &&  z.im == 0.0);

    auto w = Complex!real(1.0, 1.0);
    z = w;
    assert(z == w);
    assert(z.re == 1.0  &&  z.im == 1.0);

    auto c = Complex!float(2.0, 2.0);
    z = c;
    assert(z == c);
    assert(z.re == 2.0  &&  z.im == 2.0);
}


/*  Makes Complex!(Complex!T) fold to Complex!T.

    The rationale for this is that just like the real line is a
    subspace of the complex plane, the complex plane is a subspace
    of itself.  Example of usage:
    ---
    Complex!T addI(T)(T x)
    {
        return x + Complex!T(0.0, 1.0);
    }
    ---
    The above will work if T is both real and complex.
*/
template Complex(T)
if (is(T R == Complex!R))
{
    alias Complex = T;
}

@safe pure nothrow unittest
{
    static assert(is(Complex!(Complex!real) == Complex!real));

    Complex!T addI(T)(T x)
    {
        return x + Complex!T(0.0, 1.0);
    }

    auto z1 = addI(1.0);
    assert(z1.re == 1.0 && z1.im == 1.0);

    enum one = Complex!double(1.0, 0.0);
    auto z2 = addI(one);
    assert(z1 == z2);
}

// std.math library function overloads (KD, 2018)
@nogc
bool isNaN(Complex!double z) @safe pure nothrow
{
    double x = z.re;
    if (std.math.isNaN(x))
        return true;
    else
        return false;
}
 
@nogc
Complex!double pow(Complex!double z, int w) @safe pure nothrow
{
    Complex!double p  = complex(z.re, z.im);
    foreach ( i; 0..w-1) z *= p;
    
    return z; 
}
        
@nogc
Complex!double pow(Complex!double z, Complex!double w) @safe pure nothrow
{
    double a = z.re; double b = z.im;
    double c = w.re; double d = w.im;
    double r = std.math.sqrt(a*a + b*b);
    double theta = arg(z); Complex!double i = complex(0.0, 1.0);
    double logr = std.math.log(to!double(r));
    return exp(logr*w+i*theta*w); 
}

@nogc
Complex!double pow(Complex!double z, double w) @safe pure nothrow
{
    double a = z.re; double b = z.im;
    double c = w; double d = 0.0;
    double r = std.math.sqrt(a*a + b*b);
    double theta = arg(z); Complex!double i = complex(0.0, 1.0);
    double logr = std.math.log(to!double(r));
    return exp(logr*w+i*theta*w);
}

@nogc
Complex!double pow(double z, Complex!double w) @safe pure nothrow 
{
    double a = z; double b = 0.0;
    double c = w.re; double d = w.im;
    double r = std.math.sqrt(a*a + b*b);
    double theta = arg(complex(a, b)); Complex!double i = complex(0.0, 1.0);
    double logr = std.math.log(to!double(r));
    return exp(logr*w+i*theta*w);
}

@nogc
Complex!double fabs(Complex!double z) @safe pure nothrow 
{
    // The standard library abs() function does not satisfy analyticity, hence will not yield correct sensitivity
    // information when used in the flow solver. Below is an implementation that imposes analyticity, referenced from
    // An Automated Method for Sensitivity Analysis using Complex Variables (Martins et al, 2000).
    // A thorough explanation of the reasoning behind this implementation is provided in Martins' thesis,
    // A Coupled-Adjoint Method For High-Fidelity Aero-Structural Optimization (pg. 42, 2003).
    
    double x = z.re;
    if ( x < 0.0)
        return -z;
    else
        return z; 
}

@nogc
Complex!double exp(Complex!double z) @safe pure nothrow 
{
    double x = z.re; double y = z.im;
    double e = std.math.exp(x);
    return e*complex(std.math.cos(y), std.math.sin(y));
}

@nogc
Complex!double sqrt(Complex!double z) @safe pure nothrow 
{
    double x = z.re; double y = z.im;
    double zarg = arg(z);
    double zabs = std.math.sqrt(x*x + y*y);
    return std.math.sqrt(zabs)*complex( std.math.cos(zarg/2.0), std.math.sin(zarg/2.0) );
}

@nogc
Complex!double sin(Complex!double z) @safe pure nothrow 
{
    // The definition provided in ref.:
    // An Automated Method for Sensitivity Analysis using Complex Variables (Martins et al, 2000)
    // exhibits some error when evaluating the second example derivative computation from ref.:
    // Using Complex Variables to Estimate Derivatives of Real Functions (Squire & Trapp, 1998)
    // despite passing it's own unittest.
    // Hence we use an alternate ref.: https://proofwiki.org/wiki/Sine_of_Complex_Number
    Complex!double i = complex(0, 1);
    double a = z.re; double b = z.im;
    return std.math.sin(a) * std.math.cosh(b) + i*std.math.cos(a)*std.math.sinh(b);
}

@nogc
Complex!double cos(Complex!double z) @safe pure nothrow 
{
    // Use same ref. as sin() for consistency: https://proofwiki.org/wiki/Cosine_of_Complex_Number
    Complex!double i = complex(0, 1);
    double a = z.re; double b = z.im;
    return std.math.cos(a) * std.math.cosh(b) - i*std.math.sin(a)*std.math.sinh(b);
 }

@nogc
Complex!double tan(Complex!double z) @safe pure nothrow 
{
    double x = z.re; double y = z.im;
    Complex!double numer;
    numer = complex( std.math.sin(x)*std.math.cosh(y), std.math.cos(x)*std.math.sinh(y));
    Complex!double denom;
    denom = complex( std.math.cos(x)*std.math.cosh(y), -std.math.sin(x)*std.math.sinh(y));
    return numer/denom; 
}

@nogc
Complex!double log(Complex!double z) @safe pure nothrow 
{
    double x = z.re; double y = z.im;
    double zabs = std.math.sqrt(x*x + y*y);
    return complex( to!double(std.math.log(zabs)), arg(z));
}

@nogc
Complex!double log10(Complex!double z) @safe pure nothrow 
{
    double ln10 = std.math.log(10.0);
    return log(z)/ln10;
}

@nogc
Complex!double sinh(Complex!double z) @safe pure nothrow
{
    double x = z.re; double y = z.im;
    return (exp(z) - exp(-z))/2.0;
}

@nogc
Complex!double cosh(Complex!double z) @safe pure nothrow
{
    double x = z.re; double y = z.im;
    return (exp(z) + exp(-z))/2.0;
}

@nogc
Complex!double tanh(Complex!double z) @safe pure nothrow 
{
    Complex!double zsinh = sinh(z);
    Complex!double zcosh = cosh(z);
    return zsinh/zcosh;
}

@nogc
Complex!double fmax(Complex!double z1, Complex!double z2) @safe pure nothrow
{
        double x1 = z1.re; double x2 = z2.re;
    if (x1 >= x2)
        return z1;
    else
        return z2;
}

@nogc
Complex!double fmax(Complex!double z1, double z2) @safe pure nothrow
{
    double x1 = z1.re; double x2 = z2;
    if (x1 >= x2)
        return z1;
    else
        return complex(z2);
}

@nogc
Complex!double fmax(double z1, Complex!double z2) @safe pure nothrow 
{
    double x1 = z1; double x2 = z2.re;
    if (x1 >= x2)
        return complex(z1);
    else
        return z2;
}


@nogc
Complex!double fmin(Complex!double z1, Complex!double z2) @safe pure nothrow 
{
    double x1 = z1.re; double x2 = z2.re;
    if (x1 <= x2)
        return z1;
    else
        return z2;
}

@nogc
Complex!double fmin(double z1, Complex!double z2) { @safe pure nothrow 
    double x1 = z1; double x2 = z2.re;
    if (x1 <= x2)
        return complex(z1);
    else
        return z2;
}

@nogc
Complex!double fmin(Complex!double z1, double z2) @safe pure nothrow  
{
    double x1 = z1.re; double x2 = z2;
    if (x1 <= x2)
        return z1;
    else
        return complex(z2);
}

/*
@nogc
Complex!double sgn(Complex!double z1, Complex!double z2) @safe pure nothrow 
{
    double x1 = z1.re; double x2 = z2.re;
    if (x2 >= 0.0)
        return complex( +std.math.abs(z1.re), z1.im);
    else
        return complex( -std.math.abs(z1.re), z1.im);
}
*/

@nogc
Complex!double copysign(Complex!double z1, Complex!double z2) @safe pure nothrow 
{
    double x1 = z1.re; double x2 = z2.re;
    if (x2 >= 0.0)
        return complex( +std.math.abs(z1.re), z1.im);
    else
        return complex( -std.math.abs(z1.re), z1.im);
}

@nogc
Complex!double copysign(double z1, Complex!double z2) @safe pure nothrow 
{
    double x1 = z1.re; double x2 = z2.re;
    if (x2 >= 0.0)
        return complex( +std.math.abs(z1), 0.0);
    else
        return complex( -std.math.abs(z1), 0.0);
}

@nogc
Complex!double asin(Complex!double z) @safe pure nothrow 
{
    Complex!double i = complex(0.0, 1.0);
    return -i * log( i*z + sqrt(1.0-z*z) );
}

@nogc
Complex!double acos(Complex!double z) @safe pure nothrow 
{
    // Corrected with +ve i
    // An Automated Method for Sensitivity Analysis using Complex Variables (Martins et al, 2000).
    Complex!double i = complex(0.0, 1.0);
    return i * log( z + sqrt(z*z-1.0) );
}

@nogc
Complex!double atan(Complex!double z) @safe pure nothrow 
{
    Complex!double i = complex(0.0, 1.0);
    return 1.0/(2.0*i) * log( (i-z) / (i+z) );
}

@nogc
Complex!double atan2(Complex!double z, Complex!double w) @safe pure nothrow 
{
    // ref.: https://www.medcalc.org/manual/atan2_function.php - extension of this method to complex numbers proves problematic.
    // Below implementation provided from WolframAlpha
    Complex!double i = complex(0.0, 1.0);
    return -i * log( ((w+i*z))/(sqrt((z*z+w*w))) );
}

// end of overloaded function additions (KD, 2018)

// removed (KD, 2018)
/**
   Params: z = A complex number.
   Returns: The absolute value (or modulus) of `z`.
*/
/*
T abs(T)(Complex!T z) @safe pure nothrow @nogc
{
    import std.math : hypot;
    return hypot(z.re, z.im);
}
*/
///
@safe pure nothrow unittest
{
    static import std.math;
    assert(abs(complex(1.0)) == 1.0);
    assert(abs(complex(0.0, 1.0)) == 1.0);
    assert(abs(complex(1.0L, -2.0L)) == std.math.sqrt(5.0L));
}


/++
   Params:
    z = A complex number.
    x = A real number.
   Returns: The squared modulus of `z`.
   For genericity, if called on a real number, returns its square.
+/
T sqAbs(T)(Complex!T z) @safe pure nothrow @nogc
{
    return z.re*z.re + z.im*z.im;
}

///
@safe pure nothrow unittest
{
    import std.math;
    assert(sqAbs(complex(0.0)) == 0.0);
    assert(sqAbs(complex(1.0)) == 1.0);
    assert(sqAbs(complex(0.0, 1.0)) == 1.0);
    assert(approxEqual(sqAbs(complex(1.0L, -2.0L)), 5.0L));
    assert(approxEqual(sqAbs(complex(-3.0L, 1.0L)), 10.0L));
    assert(approxEqual(sqAbs(complex(1.0f,-1.0f)), 2.0f));
}

/// ditto
T sqAbs(T)(T x) @safe pure nothrow @nogc
if (isFloatingPoint!T)
{
    return x*x;
}

@safe pure nothrow unittest
{
    import std.math;
    assert(sqAbs(0.0) == 0.0);
    assert(sqAbs(-1.0) == 1.0);
    assert(approxEqual(sqAbs(-3.0L), 9.0L));
    assert(approxEqual(sqAbs(-5.0f), 25.0f));
}


/**
 Params: z = A complex number.
 Returns: The argument (or phase) of `z`.
 */
T arg(T)(Complex!T z) @safe pure nothrow @nogc
{
    import std.math : atan2;
    return atan2(z.im, z.re);
}

///
@safe pure nothrow unittest
{
    import std.math;
    assert(arg(complex(1.0)) == 0.0);
    assert(arg(complex(0.0L, 1.0L)) == PI_2);
    assert(arg(complex(1.0L, 1.0L)) == PI_4);
}


/**
  Params: z = A complex number.
  Returns: The complex conjugate of `z`.
*/
Complex!T conj(T)(Complex!T z) @safe pure nothrow @nogc
{
    return Complex!T(z.re, -z.im);
}

///
@safe pure nothrow unittest
{
    assert(conj(complex(1.0)) == complex(1.0));
    assert(conj(complex(1.0, 2.0)) == complex(1.0, -2.0));
}


/**
  Constructs a complex number given its absolute value and argument.
  Params:
    modulus = The modulus
    argument = The argument
  Returns: The complex number with the given modulus and argument.
*/
Complex!(CommonType!(T, U)) fromPolar(T, U)(T modulus, U argument)
    @safe pure nothrow @nogc
{
    import std.math : sin, cos;
    return Complex!(CommonType!(T,U))
        (modulus*cos(argument), modulus*sin(argument));
}

///
@safe pure nothrow unittest
{
    import std.math;
    auto z = fromPolar(std.math.sqrt(2.0), PI_4);
    assert(approxEqual(z.re, 1.0L, real.epsilon));
    assert(approxEqual(z.im, 1.0L, real.epsilon));
}

// removed (KD, 2018)
/**
    Trigonometric functions on complex numbers.

    Params: z = A complex number.
    Returns: The sine and cosine of `z`, respectively.
*/
/*
Complex!T sin(T)(Complex!T z)  @safe pure nothrow @nogc
{
    auto cs = expi(z.re);
    auto csh = coshisinh(z.im);
    return typeof(return)(cs.im * csh.re, cs.re * csh.im);
}
*/
///
/*
@safe pure nothrow unittest
{
    static import std.math;
    assert(sin(complex(0.0)) == 0.0);
    assert(sin(complex(2.0L, 0)) == std.math.sin(2.0L));
}
*/
/*
/// ditto
Complex!T cos(T)(Complex!T z)  @safe pure nothrow @nogc
{
    auto cs = expi(z.re);
    auto csh = coshisinh(z.im);
    return typeof(return)(cs.re * csh.re, - cs.im * csh.im);
}

///
@safe pure nothrow unittest
{
    import std.complex;
    assert(cos(complex(0.0)) == 1.0);
}

deprecated
@safe pure nothrow unittest
{
    import std.math;
    assert(cos(complex(0, 5.2L)) == cosh(5.2L));
    assert(cos(complex(1.3L)) == std.math.cos(1.3L));
}
*/
/**
    Params: y = A real number.
    Returns: The value of cos(y) + i sin(y).

    Note:
    `expi` is included here for convenience and for easy migration of code
    that uses $(REF _expi, std,math).  Unlike $(REF _expi, std,math), which uses the
    x87 $(I fsincos) instruction when possible, this function is no faster
    than calculating cos(y) and sin(y) separately.
*/

Complex!real expi(real y)  @trusted pure nothrow @nogc
{
    import std.math : cos, sin;
    return Complex!real(cos(y), sin(y));
}

///
@safe pure nothrow unittest
{
    import std.math : cos, sin;
    assert(expi(0.0L) == 1.0L);
    assert(expi(1.3e5L) == complex(cos(1.3e5L), sin(1.3e5L)));
}

deprecated
@safe pure nothrow unittest
{
    static import std.math;

    assert(expi(1.3e5L) == complex(std.math.cos(1.3e5L), std.math.sin(1.3e5L)));
    auto z1 = expi(1.234);
    auto z2 = std.math.expi(1.234);
    assert(z1.re == z2.re && z1.im == z2.im);
}

/**
    Params: y = A real number.
    Returns: The value of cosh(y) + i sinh(y)

    Note:
    `coshisinh` is included here for convenience and for easy migration of code
    that uses $(REF _coshisinh, std,math).
*/
Complex!real coshisinh(real y) @safe pure nothrow @nogc
{
    static import std.math;
    if (std.math.fabs(y) <= 0.5)
        return Complex!real(std.math.cosh(y), std.math.sinh(y));
    else
    {
        auto z = std.math.exp(y);
        auto zi = 0.5 / z;
        z = 0.5 * z;
        return Complex!real(z + zi, z - zi);
    }
}

///
@safe pure nothrow @nogc unittest
{
    import std.math : cosh, sinh;
    assert(coshisinh(3.0L) == complex(cosh(3.0L), sinh(3.0L)));
}

deprecated
@safe pure nothrow @nogc unittest
{
    static import std.math;
    assert(coshisinh(3.0L) == complex(std.math.cosh(3.0L), std.math.sinh(3.0L)));
    auto z1 = coshisinh(1.234);
    auto z2 = std.math.coshisinh(1.234);
    static if (real.mant_dig == 53)
    {
        assert(std.math.feqrel(z1.re, z2.re) >= real.mant_dig - 1 &&
               std.math.feqrel(z1.im, z2.im) >= real.mant_dig - 1);
    }
    else
    {
        assert(z1.re == z2.re && z1.im == z2.im);
    }
}

// removed (KD, 2018)
/**
    Params: z = A complex number.
    Returns: The square root of `z`.
*/
/*
Complex!T sqrt(T)(Complex!T z)  @safe pure nothrow @nogc
{
    static import std.math;
    typeof(return) c;
    real x,y,w,r;

    if (z == 0)
    {
        c = typeof(return)(0, 0);
    }
    else
    {
        real z_re = z.re;
        real z_im = z.im;

        x = std.math.fabs(z_re);
        y = std.math.fabs(z_im);
        if (x >= y)
        {
            r = y / x;
            w = std.math.sqrt(x)
                * std.math.sqrt(0.5 * (1 + std.math.sqrt(1 + r * r)));
        }
        else
        {
            r = x / y;
            w = std.math.sqrt(y)
                * std.math.sqrt(0.5 * (r + std.math.sqrt(1 + r * r)));
        }

        if (z_re >= 0)
        {
            c = typeof(return)(w, z_im / (w + w));
        }
        else
        {
            if (z_im < 0)
                w = -w;
            c = typeof(return)(z_im / (w + w), w);
        }
    }
    return c;
}

///
@safe pure nothrow unittest
{
    static import std.math;
    assert(sqrt(complex(0.0)) == 0.0);
    assert(sqrt(complex(1.0L, 0)) == std.math.sqrt(1.0L));
    assert(sqrt(complex(-1.0L, 0)) == complex(0, 1.0L));
}

@safe pure nothrow unittest
{
    import std.math : approxEqual;

    auto c1 = complex(1.0, 1.0);
    auto c2 = Complex!double(0.5, 2.0);

    auto c1s = sqrt(c1);
    assert(approxEqual(c1s.re, 1.09868411));
    assert(approxEqual(c1s.im, 0.45508986));

    auto c2s = sqrt(c2);
    assert(approxEqual(c2s.re, 1.1317134));
    assert(approxEqual(c2s.im, 0.8836155));
}
*/
// Issue 10881: support %f formatting of complex numbers
@safe unittest
{
    import std.format : format;

    auto x = complex(1.2, 3.4);
    assert(format("%.2f", x) == "1.20+3.40i");

    auto y = complex(1.2, -3.4);
    assert(format("%.2f", y) == "1.20-3.40i");
}

@safe unittest
{
    // Test wide string formatting
    import std.format;
    wstring wformat(T)(string format, Complex!T c)
    {
        import std.array : appender;
        auto w = appender!wstring();
        auto n = formattedWrite(w, format, c);
        return w.data;
    }

    auto x = complex(1.2, 3.4);
    assert(wformat("%.2f", x) == "1.20+3.40i"w);
}

@safe unittest
{
    // Test ease of use (vanilla toString() should be supported)
    assert(complex(1.2, 3.4).toString() == "1.2+3.4i");
}

// end of version(complex_numbers)  

} else {
    // Presume that we are building for double_numbers
    // and we do not need all the complex numbers machinery,
    // just enough to define Complex!double

import std.traits;    

struct Complex(T)
if (isFloatingPoint!T)
{
    T re;
    T im;
}

// end double_numbers version
}

version(complex_number_test) {
    import util.msg_service;
    import nm.number;
    int main() {
        // Try the various constructors and assignments.
        // PJ 2018-06-02
        Complex!double za = Complex!double("1.0");
        assert(za.re == 1.0, failedUnitTest());
        assert(za.im == 0.0, failedUnitTest());
        Complex!double zb = to!(Complex!double)("2.0");
        assert(zb.re == 2.0, failedUnitTest());
        assert(zb.im == 0.0, failedUnitTest());
        Complex!double zc = 3.0;
        assert(zc.re == 3.0, failedUnitTest());
        assert(zc.im == 0.0, failedUnitTest());
        Complex!double zd = 4;
        assert(zd.re == 4.0, failedUnitTest());
        assert(zd.im == 0.0, failedUnitTest());
        Complex!double ze = "5.0";
        assert(ze.re == 5.0, failedUnitTest());
        assert(ze.im == 0.0, failedUnitTest());
        
        // complex number reference solutions from Python 2.7.12 cmath library.
        //
        // define some test values
        Complex!double  z = complex(1.2, -3.4);
        Complex!double  w = complex(-5.3, 1.0); 
        double p = 5.1;
              
        // opCmp tests
        // Complex Cmp Complex
        assert( (z > w), failedUnitTest());
        assert( (w < z), failedUnitTest());
        assert( (z != w), failedUnitTest());

        // Complex Cmp Double
        assert( (z < p), failedUnitTest());
        assert( (p > z), failedUnitTest());
        assert( (z != p), failedUnitTest());

        // Exponent tests
        Complex!double result;
        Complex!double cpow;

        // pow(Complex, Complex)
        result = complex(0.00017039838981580696, 0.00382344206618444);
        cpow = pow(z, w);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^Complex
        cpow = z^^(w);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(Complex, Double)        
        result = complex(-10.12, -8.16);
        cpow = pow(z, 2.0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^double
        cpow = z^^2.0;
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^int
        cpow = z^^(2);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(complex, int)
        cpow = pow(z, 2);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(Double, Complex)
        result = complex(-1.6253264602682924, -1.6236827093503579);
        cpow = pow(2.0, z);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Double^^Complex
        cpow = 2.0^^z;
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Absolute value test
        Complex!double cfabs;
        // +ve real
        result = complex(1.2, -3.4);
        cfabs = fabs(z);
        assert(approxEqualNumbers(cfabs, result), failedUnitTest());

        // -ve real
        result = complex(5.3, -1.0);
        cfabs = fabs(w);
        assert(approxEqualNumbers(cfabs, result), failedUnitTest());

        // Exponential test
        result = complex(-3.209883040054176, 0.8484263372940289);
        Complex!double cexp = exp(z);
        assert(approxEqualNumbers(cexp, result), failedUnitTest());
 
        // Square root tests
        result = complex(1.5500889128472581, -1.096711282759503);
        Complex!double csqrt;
        // sqrt
        csqrt = sqrt(z);
        assert(approxEqualNumbers(csqrt, result), failedUnitTest());

        // z^^0.5
        csqrt = z^^0.5;
        assert(approxEqualNumbers(csqrt, result), failedUnitTest());
        
        // Trigonometric tests
        // sin(Complex)
        result = complex(13.979408806017995, -5.422815472463402);
        Complex!double csin = sin(z);
        assert(approxEqualNumbers(csin, result), failedUnitTest());

        // cos(Complex)
        result =complex(5.434908535625769, 13.948303613988436);
        Complex!double ccos = cos(z);
        assert(approxEqualNumbers(ccos, result), failedUnitTest());

        // tan(Complex)
        result = complex(0.0015071018758057832, -1.001642796989141);
        Complex!double ctan = tan(z);
        assert(approxEqualNumbers(ctan, result), failedUnitTest());

        // sinh(Complex)
        result = complex(-1.4593445101810318, 0.46269691906508803);
        Complex!double csinh = sinh(z);
        assert(approxEqualNumbers(csinh, result), failedUnitTest());

        // cosh(Complex)
        result = complex(-1.7505385298731442, 0.3857294182289409);
        Complex!double ccosh = cosh(z);
        assert(approxEqualNumbers(ccosh, result), failedUnitTest());
        
        // tanh(Complex)
        result = complex(0.8505969575493737, -0.0768887100657046);
        Complex!double ctanh = tanh(z);
        assert(approxEqualNumbers(ctanh, result), failedUnitTest());

        // asin(Complex)
        result = complex(0.32774305201452525, -1.990465064891069);
        Complex!double casin = asin(z);
        assert(approxEqualNumbers(casin, result), failedUnitTest());

        // acos(Complex)
        result = complex(1.2430532747803715, 1.990465064891069);
        Complex!double cacos = acos(z);
        assert(approxEqualNumbers(cacos, result), failedUnitTest());
        
        // atan(Complex)
        result = complex(1.4720985468699563, -0.2652179901713157);
        Complex!double catan = atan(z);
        assert(approxEqualNumbers(catan, result), failedUnitTest());

        // atan2(Complex, Complex)
        // cmath has no atan2 function, reference result taken from WolframAlpha
        result = complex(2.70088, 0.548252); 
        Complex!double catan2 = atan2(z, w);
        assert(approxEqualNumbers(catan2, result), failedUnitTest());
        
        // Natural log test
        result = complex(1.2824746787307684, -1.2315037123408519);
        Complex!double clog = log(z);
        assert(approxEqualNumbers(clog, result), failedUnitTest());

        // log base 10 test
        result = complex(0.5569716761534184, -0.5348352667130015);
        Complex!double clog10 = log10(z);
        assert(approxEqualNumbers(clog10, result), failedUnitTest());

        // Fmax test
        result = complex(1.2, -3.4);     
        Complex!double cmax = fmax(z, w);
        assert(approxEqualNumbers(cmax, result), failedUnitTest());

        // Fmin test
        result = complex(-5.3, 1.0);     
        Complex!double cmin = fmin(z, w);
        assert(approxEqualNumbers(cmin, result), failedUnitTest());

        // Copysign test
        result = complex(-1.2, -3.4);
        Complex!double csign = copysign(z, w);
        assert(approxEqualNumbers(csign, result), failedUnitTest());

        // Complex step finite difference derivative test
        // examples from:
        // Using Complex Variables to Estimate Derivatives of Real Functions (Squire & Trapp, 1998)

        // some common values
        Complex!double h = complex(0, 1.0e-20); // step-size
        Complex!double x0 = complex(1.5, 0.0); // original x
        Complex!double xp = x0 + h; // perturbed x

        // Example 1: F(x) = x^^(9/2) at x0 = 1.5
        double res1 = 18.600812734259759;
        Complex!double f1(Complex!double x) { return x^^(9.0/2.0); }
        double deriv1 = f1(xp).im/h.im;
        assert( std.math.approxEqual(deriv1, res1), failedUnitTest());

        // Example 2: F(x) = e^^(x)/(sin(x)^^3+cos(x)^^3) at x0 = 1.5
        double res2 = 3.62203;
        Complex!double f2(Complex!double x) { return exp(x) / ( sin(x)^^3 + cos(x)^^3 ); }
        double deriv2 = f2(xp).im/h.im;
        assert( std.math.approxEqual(deriv2, res2), failedUnitTest());
        
        return 0;
    }
}
