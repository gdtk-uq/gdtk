# zero_solvers_test.rb

$LOAD_PATH << '.'
require 'zero_solvers'
include ZeroSolvers

puts "Begin zero_solvers self-test..."
#
test_fun_1 = Proc.new do |x|
  x*x*x + x*x - 3.0*x - 3.0
end
puts ''
puts 'Test function 1.'
puts '----------------'
puts 'Example from Gerald and Wheatley, p. 45'
puts 'Solve f(x) = x^3 + x^2 - 3x -3 = 0 with initial'
puts 'guesses of x0 = 1 and x1 = 2.'
puts 'Begin function call secant()...'
puts ''
puts "Iteration \t x0 \t\tx1 \t\tx2 \t f(x2) "
puts '-----------------------------------------------------------------------'
x2 = ZeroSolvers.secant(test_fun_1, 1.0, 2.0, tol=1.0e-11,
                        limits=[], max_iterations=1000, tf=true)
puts '-----------------------------------------------------------------------'
puts 'Final result x = %g' % x2
puts 'Gerald and Wheatley report x = 1.732051'
puts 'Using bisection... x = %g' % bisection(test_fun_1, 1.0, 2.0, tol=1.0e-11)
puts ''
#
test_fun_2 = Proc.new do |x|
  3*x + Math.sin(x) - Math.exp(x)
end
puts 'Test function 2.'
puts '----------------'
puts 'Example from Gerald and Wheatley, p.45'
puts 'Solve f(x) = 3*x + sin(x) - e^x = 0 with initial'
puts 'guesses of x0 = 0 and x1 = 1.'
puts 'Begin function call secant()...'
puts ''
puts "Iteration \t x0 \t\tx1 \t\tx2 \t f(x2) "
puts '-----------------------------------------------------------------------'
x2 = ZeroSolvers.secant(test_fun_2, 0, 1, tol=1.0e-11,
                        limits=[], max_iterations=1000, tf=true)
puts '-----------------------------------------------------------------------'
puts 'Final result x = %g' % x2
puts 'Gerald and Wheatley report x = 0.3604217'
puts 'Using bisection... x = %g' % bisection(test_fun_2, 0.0, 1.0, tol=1.0e-11)
puts ''
#
test_fun_3 = Proc.new do |x|
  1.0
end
puts 'Test function 3 should throw an exception.'
begin
  x2 = ZeroSolvers.secant(test_fun_3, 0, 1)
rescue ZeroSolversError => e
  puts 'Did indeed catch the error: %s' % e.message
end
puts 'Done.'

=begin
Example transcript:

$ ruby zero_solvers_test.rb
Begin zero_solvers self-test...

Test function 1.
----------------
Example from Gerald and Wheatley, p. 45
Solve f(x) = x^3 + x^2 - 3x -3 = 0 with initial
guesses of x0 = 1 and x1 = 2.
Begin function call secant()...

Iteration 	 x0 		x1 		x2 	 f(x2) 
-----------------------------------------------------------------------
  1 	  1.000000 	 2.000000 	 1.571429 	 -1.364431e+00
  2 	  2.000000 	 1.571429 	 1.705411 	 -2.477451e-01
  3 	  1.571429 	 1.705411 	 1.735136 	 2.925540e-02
  4 	  1.705411 	 1.735136 	 1.731996 	 -5.151769e-04
  5 	  1.735136 	 1.731996 	 1.732051 	 -1.039000e-06
  6 	  1.731996 	 1.732051 	 1.732051 	 3.702993e-11
  7 	  1.732051 	 1.732051 	 1.732051 	 1.776357e-15
-----------------------------------------------------------------------
Final result x = 1.73205
Gerald and Wheatley report x = 1.732051
Using bisection... x = 1.73205

Test function 2.
----------------
Example from Gerald and Wheatley, p.45
Solve f(x) = 3*x + sin(x) - e^x = 0 with initial
guesses of x0 = 0 and x1 = 1.
Begin function call secant()...

Iteration 	 x0 		x1 		x2 	 f(x2) 
-----------------------------------------------------------------------
  1 	  1.000000 	 0.000000 	 0.470990 	 2.651588e-01
  2 	  0.000000 	 0.470990 	 0.372277 	 2.953367e-02
  3 	  0.470990 	 0.372277 	 0.359904 	 -1.294813e-03
  4 	  0.372277 	 0.359904 	 0.360424 	 5.530053e-06
  5 	  0.359904 	 0.360424 	 0.360422 	 1.021329e-09
  6 	  0.360424 	 0.360422 	 0.360422 	 -8.881784e-16
-----------------------------------------------------------------------
Final result x = 0.360422
Gerald and Wheatley report x = 0.3604217
Using bisection... x = 0.360422

Test function 3 should throw an exception.
Did indeed catch the error: Did not converge after 1000 iterations.
Done.
=end
