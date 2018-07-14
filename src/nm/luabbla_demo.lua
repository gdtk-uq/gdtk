-- luabbla_demo.lua
-- A demonstration script for using the Matrix class from
-- the bare-bones linear algebra module.
-- It can be exercised with a command like:
-- $ e4shared --custom-post --script-file=luabbla_demo.lua

print("Basic constructor resulting in nan values:")
a = Matrix:new{n=4}
assert(a:nrows() == 4)
assert(a:ncols() == 4)
print("Note that a Matrix is considered an array or rows.")
print("a=", a)

print("Simple constructor and set an element:")
a:zeros()
a:set(0,0, 5.0)
assert(a:get(0,0) == 5.0)
print("Note that indices start at zero.")
print("a=", a)

print("Some arithmetic:")
b = Matrix:new{nrows=4, ncols=4}
b:zeros()
c = b + a
assert(c:get(0,0) == 5.0)
print("c=", c)

d = b - a
assert(d:get(0,0) == -5.0)
print("d=", d)

print("Set to identity matrix:")
a:eye()
assert(a:get(1,1) == 1.0)
print("a=", a)

print("Construct a rectangular matrix:")
e = zeros(4, 2)
assert(e:get(3,1) == 0.0)
print("e=", e)

print("Identity constructor:")
f = eye(5)
assert(f:get(2,2) == 1.0)

print("Set and solve a matrix equation:")
A = Matrix:new{n=4}
A:set(0,0, 0.0); A:set(0,1, 2.0); A:set(0,2, 0.0); A:set(0,3, 1.0); 
A:set(1,0, 2.0); A:set(1,1, 2.0); A:set(1,2, 3.0); A:set(1,3, 2.0); 
A:set(2,0, 4.0); A:set(2,1, -3.0); A:set(2,2, 0.0); A:set(2,3, 1.0); 
A:set(3,0, 6.0); A:set(3,1, 1.0); A:set(3,2, -6.0); A:set(3,3, -5.0); 
b = Matrix:new{vec={0.0, -2.0, -7.0, 6.0}, orient='column'}
print("A=", A)
print("b=", b)
x = solve(A, b)
assert(math.abs(x:get(0,0) - -0.5) < 1.0e-6)
assert(math.abs(x:get(1,0) - 1.0) < 1.0e-6)
assert(math.abs(x:get(2,0) - 0.333333) < 1.0e-6)
assert(math.abs(x:get(3,0) - -2.0) < 1.0e-6)
print("A=", A)
print("b=", b)
print("x=", x)

print("Stack some matrices horizontally:")
a1 = eye(2)
a2 = zeros(2, 2)
a3 = eye(2)
h = hstack{a1, a2, a3}
print("h=", h)

print("Build augmented matrix:")
Ab = hstack{A, eye(4)}
print("Ab=", Ab)
print("Apply Gauss-Jordan elimination:")
gaussJordan(Ab)
print("Ab=", Ab)

print("Extract inverse matrix:")
Ainv = Ab:getSlice(0,4, 4,8)
print("Ainv=", Ainv)
print("Compute solution via dot product and compare:")
x2 = Ainv:dot(b)
print("x2=", x2)
diff = x2 - x
print("diff=", diff)

print("Done.")
