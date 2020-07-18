# add-time-column.awk
# Get the T4 experimental pressure recording into desired format.
# PJ, 2020-05-26
BEGIN {
    print("# t(us) p(kPa)")
    dt = 1.0 # microseconds
    t = 0.0
}
$1 != "#" && $1 != "" {
    print(t, $1)
    t = t + dt
}
