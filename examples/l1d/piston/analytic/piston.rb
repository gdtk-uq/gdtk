#!/usr/bin/env ruby
# piston.rb
# Piston motion calculation, assuming ideal conditions.
# For theory, see PJ Fluid Mechanics Lecture notes, 1994,
# and AGARDograph 138 Ballistic Range Technology, 1970.
# Dimensional quantities from Rowan's case 2 in Casbar paper, 2007.
#
# PJ, 2019-05-17
#
# Driving gas is ideal air
R_air = 287.10 # J/(kg.K)
G = 1.4
# Initial condition for gas.
P_0 = 1.0e5 # Pa
rho_0 = 1.0 # kg/m**3
A_0 = Math.sqrt(G*P_0/rho_0) # m/s
# Piston and tube dimensions.
Area = Math::PI*0.2**2 # m**2
M_p = 1.0 # kg

def t_bar(t)
  P_0*Area*t/(M_p*A_0)
end

def x(t)
  tb = t_bar(t)
  t1 = 2.0/(G+1)
  t2 = 1 + tb*(G+1)/2
  t3 = 1 + tb - t2**t1
  x_bar = 2/(G-1) * t3
  # and back to dimensional form
  x_bar*M_p*(A_0**2)/(P_0*Area)
end

def v(t)
  t1 = (1-G)/(G+1)
  t2 = 1 + t_bar(t)*(G+1)/2
  t3 = 1 - t2**t1
  v_bar = 2/(G-1) * t3
  # and back to dimensional form
  v_bar*A_0
end

puts "# Computed quantities: A_0=#{A_0}, Area=#{Area}"
puts "# t(s)   x(m)   v(m/s)"
dt = 0.001 # sec
(0..40).each do |i|
  t = i*dt
  puts "#{t} #{x(t)} #{v(t)}"
end
