BEGIN {
    print "# 1:iteration 2:wall-clock[s] 3:mass-residual 4:x-mom-residual 5:y-mom-residual 6:z-mom-residual 7:energy residual"
}
/^RESIDUALS/ {
    print $3, $5, $7, $9, $11, $13, $15
}

