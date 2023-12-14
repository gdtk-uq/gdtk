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

    --------------------------------------------------------------------------------------
    This module has been adapted for use in complexifying our multi-physics fluid solver.
    We use this module as the basis to perform complex variable differentiation to construct
    numerical Jacobians and design sensitivities for shape optimization. For more information
    on complex variable differentiation, refer to:

    [1] Squire et al., Using Complex Variables to Estimate Derivatives of Real Functions,
        SIAM review, vol. 40, No. 1, pp. 110-112, March 1998.

    I have marked all modifications made to the original module and provided an accompanying
    description or reference. For further information on developing a complex number module for
    complex variable differentation, refer to:

    [1] Martins et al., An Automated Method for Sensitivity Analysis using Complex Variables,
        38th Aerospace Sciences Meeting and Exhibit, AIAA paper 2000-0689, 2000.

    [2] Martins, A Coupled-Adjoint Method for High-Fidelity Aero-Structural Optimization,
        Stanford University, 2002.

    [3] Nielsen et al., Efficient Construction of Discrete Adjoint Operators on Unstructured Grids
        by Using Complex Variables, AIAA Journal, Vol. 44, No. 4, 2006, pp. 827–836.

    author: Kyle A. Damm
    date:   2018 (modified 2020 and 2023)
    --------------------------------------------------------------------------------------

*/
module ntypes.complex;

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

    // CASTING OPERATORS

    R opCast(R : T)() const
    {
        return re;
    }

    // COMPARISON OPERATORS

    // For the opEquals comparison operators, we only want the real component to be checked,
    // the imaginary component will most likely (and necessarily) have a different value, since
    // it encodes the sensitivity information which is propagated forward through the calculations.
    // KAD 2018

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

    // added an OpCmp comparison operator
    // implemented from sources:
    // + https://forum.dlang.org/post/p1qa9n$2hje$1@digitalmars.com
    // + http://www.angelcode.com/angelscript/sdk/docs/manual/doc_script_class_ops.html
    // KAD 2018

    // complex version
    int opCmp(Complex!double z) const
    {
        import std.math: fabs;
        auto diff = re - z.re;
        auto epsilon = 1.0e-50;
        if (  fabs(diff) < epsilon )
            return 0;
        else if ( diff < 0 )
            return -1;
        else
            return 1;
    }

    // numeric version
    int opCmp(double z) const
    {
        import std.math: fabs;
        auto diff = re - z.re;
        auto epsilon = 1.0e-50;
        if (  fabs(diff) < epsilon )
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
    ref Complex opOpAssign(string op, C)(C z)
        if (op == "^^" && is(C R == Complex!R))
    {
        import std.math : exp, log, cos, sin;
        immutable r = cabs(this);
        immutable t = arg(this);
        immutable ab = r^^z.re * exp(-t*z.im);
        immutable ar = t*z.re + log(r)*z.im;

        re = ab*cos(ar);
        im = ab*sin(ar);
        return this;
    }

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

    // complex ^^= real
    ref Complex opOpAssign(string op, R)(R r)
        if (op == "^^" && isFloatingPoint!R)
    {
        import std.math : cos, sin;
        if (isInteger(r)) {
            int p = cast(int) r;
            this ^^= p;
            return this;
        }
        immutable ab = cabs(this)^^r;
        immutable ar = arg(this)*r;
        re = ab*cos(ar);
        im = ab*sin(ar);
        return this;
    }

    // complex ^^= int
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
            // Robb Watt found that for the instances when re is some
            // small negative number and im is 0.0 calling the ^^= real
            // opOpAssign was returning a value with a non-zero im component.
            // For int powers we are able to perform an explicit loop, which
            // appears to be more robust. See pow(Complex!T, int). KAD 2023.
            Complex!T p = complex(this.re, this.im);
            foreach (k; 0..abs(i)-1) this *= p;
            if (i < 0) this = 1.0/this;
        }
        return this;
    }

}

@("Test casting")
@safe pure nothrow unittest
{
    auto c1 = complex(1.0, 1.0);
    auto c2 = cast(double) c1;
    assert(is(typeof(c2) == double));
    assert(c2 == 1.0);
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
    assert(isClose(cabs(ctc), cabs(c1)*cabs(c2), EPS));
    assert(isClose(arg(ctc), arg(c1)+arg(c2), EPS));

    auto cdc = c1 / c2;
    assert(isClose(cabs(cdc), cabs(c1)/cabs(c2), EPS));
    assert(isClose(arg(cdc), arg(c1)-arg(c2), EPS));

    auto cec = c1^^c2;
    assert(isClose(cec.re, 0.11524131979943839881, EPS));
    assert(isClose(cec.im, 0.21870790452746026696, EPS));

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
    assert(isClose(cabs(cdr), cabs(c1)/a, EPS));
    assert(isClose(arg(cdr), arg(c1), EPS));

    auto cer = c1^^3.0;
    assert(isClose(cabs(cer), cabs(c1)^^3, EPS));
    assert(isClose(arg(cer), arg(c1)*3, EPS));

    auto rpc = a + c1;
    assert(rpc == cpr);

    auto rmc = a - c1;
    assert(rmc.re == a-c1.re);
    assert(rmc.im == -c1.im);

    auto rtc = a * c1;
    assert(rtc == ctr);

    auto rdc = a / c1;
    assert(isClose(cabs(rdc), a/cabs(c1), EPS));
    assert(isClose(arg(rdc), -arg(c1), EPS));

    rdc = a / c2;
    assert(isClose(cabs(rdc), a/cabs(c2), EPS));
    assert(isClose(arg(rdc), -arg(c2), EPS));

    auto rec1a = 1.0 ^^ c1;
    assert(rec1a.re == 1.0);
    assert(rec1a.im == 0.0);

    auto rec2a = 1.0 ^^ c2;
    assert(rec2a.re == 1.0);
    assert(rec2a.im == 0.0);

    auto rec1b = (-1.0) ^^ c1;
    assert(isClose(cabs(rec1b), std.math.exp(-PI * c1.im), EPS));
    auto arg1b = arg(rec1b);
    /* The argument _should_ be PI, but floating-point rounding error
     * means that in fact the imaginary part is very slightly negative.
     */
    assert(isClose(arg1b, PI, EPS) || isClose(arg1b, -PI, EPS));

    auto rec2b = (-1.0) ^^ c2;
    assert(isClose(cabs(rec2b), std.math.exp(-2 * PI), EPS));
    assert(isClose(arg(rec2b), PI_2, EPS));

    auto rec3a = 0.79 ^^ complex(6.8, 5.7);
    auto rec3b = complex(0.79, 0.0) ^^ complex(6.8, 5.7);
    assert(isClose(rec3a.re, rec3b.re, EPS));
    assert(isClose(rec3a.im, rec3b.im, EPS));

    auto rec4a = (-0.79) ^^ complex(6.8, 5.7);
    auto rec4b = complex(-0.79, 0.0) ^^ complex(6.8, 5.7);
    assert(isClose(rec4a.re, rec4b.re, 0.0, EPS));
    assert(isClose(rec4a.im, rec4b.im, 0.0, EPS));

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
    // NOTE: The arg is approx -2.44e-16, EPS = 2.22e-16
    assert(isClose(arg(rer2), 0.0, 0.0, 2*EPS));

    auto rer3 = (-a) ^^ complex(-2.0, 0.0);
    rcheck = (-a) ^^ (-2.0);
    assert(feqrel(rer3.re, rcheck) == double.mant_dig);
    assert(isIdentical(rer3.re, rcheck));
    assert(isClose(rer3.im, 0.0, 0.0, EPS));

    auto rer4 = a ^^ complex(-2.0, 0.0);
    rcheck = a ^^ (-2.0);
    assert(feqrel(rer4.re, rcheck) == double.mant_dig);
    assert(isIdentical(rer4.re, rcheck));
    assert(rer4.im == 0.0);

    // Check Complex-int operations.
    foreach (i; 0 .. 6)
    {
        auto cei = c1^^i;
        assert(isClose(cabs(cei), cabs(c1)^^i, EPS * (i == 0 ? 1 : i)));
        // Use cos() here to deal with arguments that go outside
        // the (-pi,pi] interval (only an issue for i>3).
        assert(isClose(std.math.cos(arg(cei)), std.math.cos(arg(c1)*i), EPS * (i == 0 ? 1 : i)));
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
    assert(isClose(c1c.re, 1.0, EPS));
    assert(isClose(c1c.im, 0.0, EPS));

    c1c = c1;
    c1c /= c2;
    // WARN: The expected value precision is too low to use EPS as the comparison, 
    //       as the value is not a round number. More digits shown below 
    // assert(isClose(c1c.re, 0.58823529411764705, EPS));
    // assert(isClose(c1c.im, -0.35294117647058823, EPS));
    assert(isClose(c1c.re, 0.588235, 1E-06));
    assert(isClose(c1c.im, -0.352941, 1E-06));

    c2c /= c1;
    assert(isClose(c2c.re, 1.25, EPS));
    assert(isClose(c2c.im, 0.75, EPS));

    c2c = c2;
    c2c /= c2;
    assert(isClose(c2c.re, 1.0, EPS));
    assert(isClose(c2c.im, 0.0, EPS));
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

// Added several std.math library function overloads...
// most of these are taken from the references in the header of this file.
// The original complex module definitions for sin, cos, sqrt and abs are unmodified.
// KAD 2018, 2020, 2023
@nogc
bool isInteger(double n) @safe pure nothrow
{
    import std.math: fabs;

    // we will return false immediately if n is
    // not representable as an integer since we use
    // this routine to determine whether it is safe
    // to cast a double as an integer
    if (fabs(n) > int.max) { return false; }

    // convert float to integer
    int x = cast(int) n;

    double diff = n - x;

    // if n is not equivalent to any integer
    if (fabs(diff) > 0) { return false; }

    // else n is an integer
    return true;
}

@nogc
bool isNaN(T)(Complex!T z) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    import std.math: isNaN;
    immutable x = z.re;
    if (isNaN(x))
        return true;
    else
        return false;
}

@nogc
Complex!T pow(T)(Complex!T z, int w) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    switch (w)
    {
    case 0:
        z.re = 1.0;
        z.im = 0.0;
        break;
    case 1:
        // identity; do nothing
        break;
    case 2:
        z *= z;
        break;
    case 3:
        auto p = z;
        z *= p;
        z *= p;
        break;
    default:
        // Robb Watt found that for the instances when re is some
        // small negative number and im is 0.0 calling the ^^= real
        // opOpAssign was returning a value with a non-zero im component.
        // For int powers we are able to perform an explicit loop, which
        // appears to be more robust. See ^^= int. KAD 2023.
        Complex!T p = complex(z.re, z.im);
        foreach (k; 0..abs(w)-1) z *= p;
        if (w < 0) z = 1/z;
    }
    return z;
}

@nogc
Complex!T pow(T)(Complex!T z, Complex!T w) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    import std.math: sqrt, log;
    immutable r = sqrt(z.re*z.re + z.im*z.im);
    immutable theta = arg(z);
    Complex!double i = complex(0.0, 1.0);
    immutable logr = log(r);
    return typeof(return)(exp(logr*w+i*theta*w));
}

@nogc
Complex!T pow(T)(Complex!T z, T w) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    import std.math: sqrt, log;
    if (isInteger(w)) {
        int p = cast(int) w;
        return pow(z, p);
    }
    immutable r = sqrt(z.re*z.re + z.im*z.im);
    immutable theta = arg(z);
    Complex!double i = complex(0.0, 1.0);
    immutable logr = log(r);
    return typeof(return)(exp(logr*w+i*theta*w));
}

@nogc
Complex!T pow(T)(T z, Complex!T w) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    import std.math: sqrt, log;
    immutable theta = arg(complex(z, 0.0));
    Complex!double i = complex(0.0, 1.0);
    immutable logr = log(z);
    return typeof(return)(exp(logr*w+i*theta*w));
}

@nogc
Complex!T fabs(T)(Complex!T z) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    // A standard library implementation of the abs() function does not satisfy analyticity.
    // Below is an implementation that imposes analyticity, taken from the references in the header.
    // KAD 2018
    immutable x = z.re;
    if ( x < 0.0)
        return -z;
    else
        return z;
}

Complex!T exp(T)(Complex!T z) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    import std.math: exp, cos, sin;
    immutable e = exp(z.re);
    return e*complex(cos(z.im), sin(z.im));
}

@nogc
Complex!T tan(T)(Complex!T z) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    import std.math: cos, sin, cosh, sinh;
    immutable re = z.re;
    immutable im = z.im;
    Complex!T numer = complex( sin(re)*cosh(im), cos(re)*sinh(im));
    Complex!T denom = complex( cos(re)*cosh(im), -sin(re)*sinh(im));
    return numer/denom;
}

@nogc
Complex!T log(T)(Complex!T z) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    import std.math: sqrt, log;
    immutable re = z.re;
    immutable im = z.im;
    immutable zabs = sqrt(re*re + im*im);
    return typeof(return)(complex( log(zabs), arg(z) ));
}

@nogc
Complex!T log10(T)(Complex!T z) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    import std.math;
    immutable ln10 = std.math.log(10.0);
    return typeof(return)(log(z)/ln10);
}

@nogc
Complex!T sinh(T)(Complex!T z) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    return (exp(z) - exp(-z))/2.0;
}

@nogc
Complex!T cosh(T)(Complex!T z) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    return (exp(z) + exp(-z))/2.0;
}

@nogc
Complex!T tanh(T)(Complex!T z) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    return sinh(z)/cosh(z);
}

@nogc
Complex!T fmax(T)(Complex!T z1, Complex!T z2) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    immutable x1 = z1.re;
    immutable x2 = z2.re;
    if (x1 >= x2)
        return z1;
    else
        return z2;
}

@nogc
Complex!T fmax(T)(Complex!T z1, T z2) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    immutable x1 = z1.re;
    immutable x2 = z2;
    if (x1 >= x2)
        return z1;
    else
        return complex(z2);
}

@nogc
Complex!T fmax(T)(T z1, Complex!T z2) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    immutable x1 = z1;
    immutable x2 = z2.re;
    if (x1 >= x2)
        return complex(z1);
    else
        return z2;
}

@nogc
Complex!T fmin(T)(Complex!T z1, Complex!T z2) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    immutable x1 = z1.re;
    immutable x2 = z2.re;
    if (x1 <= x2)
        return z1;
    else
        return z2;
}

@nogc
Complex!T fmin(T)(T z1, Complex!T z2)  @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    immutable x1 = z1;
    immutable x2 = z2.re;
    if (x1 <= x2)
        return complex(z1);
    else
        return z2;
}

@nogc
Complex!T fmin(T)(Complex!T z1, T z2) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    immutable x1 = z1.re;
    immutable x2 = z2;
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
Complex!T copysign(T)(Complex!T z1, Complex!T z2) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    import std.math: abs;
    immutable x1 = z1.re;
    immutable x2 = z2.re;
    if (x2 >= 0.0)
        return complex( +abs(z1.re), z1.im);
    else
        return complex( -abs(z1.re), z1.im);
}

@nogc
Complex!T copysign(T)(T z1, Complex!T z2) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    import std.math: abs;
    immutable x1 = z1.re;
    immutable x2 = z2.re;
    if (x2 >= 0.0)
        return complex( +abs(z1), 0.0);
    else
        return complex( -abs(z1), 0.0);
}

@nogc
Complex!T asin(T)(Complex!T z) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    Complex!double i = complex(0.0, 1.0);
    return -i * log( i*z + sqrt(1.0-z*z) );
}

@nogc
Complex!T acos(T)(Complex!T z) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    // Corrected with +ve i
    // An Automated Method for Sensitivity Analysis using Complex Variables (Martins et al, 2000).
    Complex!double i = complex(0.0, 1.0);
    return i * log( z + sqrt(z*z-1.0) );
}

@nogc
Complex!T atan(T)(Complex!T z) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    Complex!double i = complex(0.0, 1.0);
    return 1.0/(2.0*i) * log( (i-z) / (i+z) );
}

@nogc
Complex!T atan2(T)(Complex!T z, Complex!T w) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    // ref.: https://www.medcalc.org/manual/atan2_function.php - extension of this method to complex numbers proves problematic.
    // Below implementation provided from WolframAlpha
    Complex!double i = complex(0.0, 1.0);
    return -i * log( ((w+i*z))/(sqrt((z*z+w*w))) );
}

// end of overloaded function additions
// KAD 2018, 2020, 2023

T cabs(T)(Complex!T z) @safe pure nothrow @nogc
{
    import std.math : hypot;
    return hypot(z.re, z.im);
}

///
@safe pure nothrow unittest
{
    static import std.math;
    assert(cabs(complex(1.0)) == 1.0);
    assert(cabs(complex(0.0, 1.0)) == 1.0);
    assert(cabs(complex(1.0L, -2.0L)) == std.math.sqrt(5.0L));
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
    assert(isClose(sqAbs(complex(1.0L, -2.0L)), 5.0L));
    assert(isClose(sqAbs(complex(-3.0L, 1.0L)), 10.0L));
    assert(isClose(sqAbs(complex(1.0f,-1.0f)), 2.0f));
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
    assert(isClose(sqAbs(-3.0L), 9.0L));
    assert(isClose(sqAbs(-5.0f), 25.0f));
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
    auto z = fromPolar(std.math.sqrt(2.0L), PI_4);
    assert(isClose(z.re, 1.0L, real.epsilon));
    assert(isClose(z.im, 1.0L, real.epsilon));
}

/**
    Trigonometric functions on complex numbers.

    Params: z = A complex number.
    Returns: The sine and cosine of `z`, respectively.
*/
Complex!T sin(T)(Complex!T z)  @safe pure nothrow @nogc
{
    auto cs = expi(z.re);
    auto csh = coshisinh(z.im);
    return typeof(return)(cs.im * csh.re, cs.re * csh.im);
}
///
@safe pure nothrow unittest
{
    static import std.math;
    assert(sin(complex(0.0)) == 0.0);
    assert(sin(complex(2.0L, 0)) == std.math.sin(2.0L));
}
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
    assert(cos(complex(0, 5.2L)) == std.math.cosh(5.2L));
    assert(cos(complex(1.3L)) == std.math.cos(1.3L));
}
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

/**
    Params: z = A complex number.
    Returns: The square root of `z`.
*/
@nogc
Complex!T sqrt(T)(Complex!T z) @safe pure nothrow
if ( is(typeof(T(0.0)) == double) ||
     is(typeof(T(0.0)) == float)  ||
     is(typeof(T(0.0)) == real))
{
    import std.math : sqrt, cos, sin;
    T x = z.re; T y = z.im;
    T zarg = arg(z);
    T zabs = std.math.sqrt(x*x + y*y);
    return sqrt(zabs)*complex( cos(zarg/2.0), sin(zarg/2.0) );
}
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
*/
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
    assert(isClose(c1s.re, 1.09868411, 1E-08));
    assert(isClose(c1s.im, 0.45508986, 1E-08));

    auto c2s = sqrt(c2);
    assert(isClose(c2s.re, 1.1317139, 1E-07));
    assert(isClose(c2s.im, 0.8836155, 1E-07));
}

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

        // test isInteger() routine
        double val1 = 10.234;
        bool result1 = isInteger(val1);
        assert( result1 == false, failedUnitTest());
        double val2 = -10.234;
        bool result2 = isInteger(val2);
        assert( result2 == false, failedUnitTest());
        double val3 = 10.0;
        bool result3 = isInteger(val3);
        assert( result3 == true, failedUnitTest());
        double val4 = -10.0;
        bool result4 = isInteger(val4);
        assert( result4 == true, failedUnitTest());
        double val5 = 10.0 + 1.0e-15;
        bool result5 = isInteger(val5);
        assert( result5 == false, failedUnitTest());
        double val6 = 10.0 + 1.0e-16;
        bool result6 = isInteger(val6);
        assert( result6 == true, failedUnitTest());
        double val7 = 2.0 * int.max;
        bool result7 = isInteger(val7);
        assert( result7 == false, failedUnitTest());

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
        result = complex(1.0, 0.0);
        cpow = pow(z, 0.0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^double
        cpow = z^^0.0;
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^int
        cpow = z^^(0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(complex, int)
        cpow = pow(z, 0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(Complex, Double)
        result = z;
        cpow = pow(z, 1.0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^double
        cpow = z^^1.0;
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^int
        cpow = z^^(1);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(complex, int)
        cpow = pow(z, 1);
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

        // pow(Complex, Double)
        result = complex(-39.888, 24.616);
        cpow = pow(z, 3.0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^double
        cpow = z^^3.0;
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^int
        cpow = z^^(3);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(complex, int)
        cpow = pow(z, 3);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(Complex, Double)
        result = complex(985.105088, -1963.766016);
        cpow = pow(z, 6.0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^double
        cpow = z^^6.0;
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^int
        cpow = z^^(6);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(complex, int)
        cpow = pow(z, 6);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(Complex, Double)
        result = complex(1.0, 0.0);
        cpow = pow(z, -0.0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^double
        cpow = z^^(-0.0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^int
        cpow = z^^(-0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(complex, int)
        cpow = pow(z, -0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(Complex, Double)
        result = complex(0.09230769230769233, 0.26153846153846155);
        cpow = pow(z, -1.0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^double
        cpow = z^^(-1.0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^int
        cpow = z^^(-1);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(complex, int)
        cpow = pow(z, -1);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(Complex, Double)
        result = complex(-0.05988165680473373, 0.04828402366863906);
        cpow = pow(z, -2.0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^double
        cpow = z^^(-2.0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^int
        cpow = z^^(-2);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(complex, int)
        cpow = pow(z, -2);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(Complex, Double)
        result = complex(-0.018155666818388715, -0.011204369594902138);
        cpow = pow(z, -3.0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^double
        cpow = z^^(-3.0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^int
        cpow = z^^(-3);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(complex, int)
        cpow = pow(z, -3);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(Complex, Double)
        result = complex(0.00020409033960117344, 0.0004068456025502563);
        cpow = pow(z, -6.0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^double
        cpow = z^^(-6.0);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // Complex^^int
        cpow = z^^(-6);
        assert(approxEqualNumbers(cpow, result), failedUnitTest());

        // pow(complex, int)
        cpow = pow(z, -6);
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
        assert( std.math.isClose(deriv1, res1), failedUnitTest());

        // Example 2: F(x) = e^^(x)/(sin(x)^^3+cos(x)^^3) at x0 = 1.5
        double res2 = 3.6220337007163259;
        Complex!double f2(Complex!double x) { return exp(x) / ( sin(x)^^3 + cos(x)^^3 ); }
        double deriv2 = f2(xp).im/h.im;
        // debug { import std.stdio; writefln("deriv2=%.18f res2=%.18f", deriv2, res2); }
        assert( std.math.isClose(deriv2, res2), failedUnitTest());

        return 0;
    }
}
