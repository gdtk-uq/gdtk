# locate_shock.awk

BEGIN {
    p_old = 0.0;
    x_old = -2.0;   # dummy position
    y_old = -2.0;
    p_trigger = 2.0e6; # something midway between free stream and stagnation
    shock_found = 0;
}

$1 != "#" { # for any non-comment line, do something
    p_new = $9;
    x_new = $1;
    y_new = $2;
    # print "p_new=", p_new, "x_new", x_new, "y_new", y_new
    if ( p_new > p_trigger && shock_found == 0 ) {
        shock_found = 1;
        frac = (p_new - p_trigger) / (p_new - p_old);
	x = x_old + frac * (x_new - x_old);
	y = y_old + frac * (y_new - y_old);
	print "shock-location= ", x, y
    }
    p_old = p_new;
    x_old = x_new;
    y_old = y_new;
}

END {
    if ( shock_found == 0 ) {
        print "shock not located";
    }
    print "done."
}
