// limiters.d

module limiters;

import std.math;
import nm.complex;
import nm.number;

@nogc void min_mod_limit(ref number a, ref number b)
// A rough slope limiter.
{
    if (a * b < 0.0) {
        a = 0.0;
        b = 0.0;
    } else {
        a = copysign(fmin(fabs(a), fabs(b)), a);
        b = a;
    }
}

@nogc void van_albada_limit(ref number a, ref number b)
// A smooth slope limiter.
{
    immutable double eps = 1.0e-12;
    number s = (a*b + fabs(a*b))/(a*a + b*b + eps);
    a *= s;
    b *= s;
}

@nogc number clip_to_limits(number q, number A, number B)
// Returns q if q is between the values A and B, else
// it returns the closer limit of the range [min(A,B), max(A,B)].
{
    number lower_limit = (A <= B) ? A : B;
    number upper_limit = (A > B) ? A : B;
    number qclipped = (q > lower_limit) ? q : lower_limit;
    return (qclipped <= upper_limit) ? qclipped : upper_limit;
}
