# show-residuals.gnuplot
#
# Run with the command
# $ gnuplot show-residual.gnuplot
#
# The array of plots will refresh periodically.
# Press 'q' or 'x' to exit.
#
# PJ, 2019-10-31
#
job = "swbli"
period_seconds = 2.0

residuals_file = "./config/" . job ."-residuals.txt"
window_title = "Eilmer4 residuals for ".job
set term qt title window_title
# set xlabel "step"
# set ylabel "residual"
set logscale y
bind "x" "end=1"
bind "q" "end=1"
end = 0
while( end==0 ) {
  set multiplot layout 2,2
  plot residuals_file using "step":"x-mom" with lines, \
    residuals_file using "step":"y-mom" with lines
  plot residuals_file using "step":"energy" with lines
  plot residuals_file using "step":"mass" with lines
  plot residuals_file using "step":"L2" with lines
  pause period_seconds
}

