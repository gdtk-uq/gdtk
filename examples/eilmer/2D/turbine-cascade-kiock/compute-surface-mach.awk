# compute-surface-mach.awk
BEGIN {
    p0 = 100.0e3;
    g = 1.4;
    gm1 = g - 1.0;
    angle = 33.3*3.14159/180.0;
    print("# chord Mach p");
}
!/pos.x/ {
    x = $1;
    y = $2;
    p = $9;
    if (p > p0) { p = p0; }; # limit to stagnation pressure
    # Compute chord position for surface point.
    c = x*cos(angle) - y*sin(angle);
    # Compute isentropic Mach number for teh surface pressure.
    Mach = sqrt(2.0/gm1*((p0/p)^(gm1/g) - 1.0));
    print(c, Mach, p);
}
