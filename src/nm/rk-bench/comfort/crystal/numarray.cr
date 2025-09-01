# numarray.cr
#
# Element-by-element vector arithmetic, like numpy.
#
# PJ 2025-08-20

require "math"

class NArray(T)
  property data : Array(T)

  def initialize(data : Array(T))
    # Note that we just keep a reference to the incoming array.
    @data = data
  end

  def self.new(n : Int, value : T)
    data = Array.new(n, value)
    NArray.new(data)
  end

  def self.new(other : NArray(T))
    NArray.new(other.data.clone)
  end

  def clone() : NArray(T)
    NArray.new(@data.clone)
  end

  def to_s(io)
    io << "["
    n = @data.size
    n.times do |i|
      io << @data[i].to_s
      if i < n-1
        io << ", "
      end
    end
    io << "]"
  end

  def [](i : Int)
    @data[i]
  end

  def []=(i : Int, value : T)
    @data[i] = value
  end

  def size()
    @data.size
  end

  def +(other : NArray(T)) : NArray(T)
    n = self.size
    if other.size != n
      raise "Mismatch in NArray sizes: self.size=#{n} other.size=#{other.size}."
    end
    result = NArray.new(self)
    n.times do |i|
      result[i] = self[i] + other[i]
    end
    result
  end

  def -(other : NArray(T)) : NArray(T)
    n = self.size
    if other.size != n
      raise "Mismatch in NArray sizes: self.size=#{n} other.size=#{other.size}."
    end
    result = NArray.new(self)
    n.times do |i|
      result[i] = self[i] - other[i]
    end
    result
  end

  def -() : NArray(T)
    n = self.size
    result = NArray.new(self)
    n.times do |i|
      result[i] = -self[i]
    end
    result
  end

  def *(other : Number) : NArray(T)
    n = self.size
    result = NArray.new(self)
    n.times do |i|
      result[i] *= other
    end
    result
  end

  def /(other : Number) : NArray(T)
    n = self.size
    result = NArray.new(self)
    n.times do |i|
      result[i] /= other
    end
    result
  end

  def copy!(other : NArray(T))
    n = self.size
    n.times do |i|
      self[i] = other[i]
    end
  end

  def abs!()
    n = self.size
    n.times do |i|
      self[i] = self[i].abs()
    end
  end

  def dot(other : NArray(T)) : T
    sum : T = 0
    n = self.size
    if other.size != n
      raise "Mismatch in NArray sizes: self.size=#{n} other.size=#{other.size}."
    end
    n.times { |i| sum += self[i] * other[i] }
    sum
  end

end # class

# We want to extend the Number class/struct so that we can multiply scalars (on the left)
# by array elements (on the right) and end up with an array result that is scaled.

struct Number

  def *(other : NArray(T)) : NArray(Float64)
    n = other.size
    result = NArray.new(other)
    n.times do |i|
      result[i] *= self
    end
    result
  end

end
