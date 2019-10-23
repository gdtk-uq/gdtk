#!/usr/bin/env ruby
# xtdata.rb
# Accumulate x,t data for the selected variable, over the tindx saved solutions.
# PJ, 2019-10-22 initial experiment for Hans H. and David G.
#     2019-10-23 Make the script a bit more versatile.
#
require 'getoptlong'
require 'open3'
require 'fileutils'

opts = GetoptLong.new(
  ["--job", "-j", GetoptLong::REQUIRED_ARGUMENT],
  ["--vcolumn", "-v", GetoptLong::REQUIRED_ARGUMENT],
  ["--xcolumn", "-x", GetoptLong::REQUIRED_ARGUMENT],
  ["--log10", "-l", GetoptLong::NO_ARGUMENT],
  ["--output", "-o", GetoptLong::REQUIRED_ARGUMENT],
  ["--slice", "-s", GetoptLong::REQUIRED_ARGUMENT],
  ["--help", "-h", GetoptLong::NO_ARGUMENT]
)

puts "xtdata"
puts "Accumulate x,t,var data over the saved solutions."

job_name = ""
var_name = ""
# Default columns select data for x and p.
x_column = 1
var_column = 9
take_log10 = false
output_name = "xt.data"
slice_str = ":,:,0,0"

opts.each do |opt, arg|
  case opt
  when /-h/
    puts "Usage:"
    puts "$ xtdata.rb --job=jobName ?--vcolumn=9? ?--xcolumn=1?"
    puts "    ?--log10? ?--output=xt.data? ?--slice=:,:,0,0?"
    puts "Remember that column counting starts at 1 in the labels for variables."
    puts "The default slice string assumes that all blocks are in a single line"
    puts "along the x-axis and, for each block, you want to select the full i-range"
    puts "of cells with j=0 and k=0 cell indices."
    exit 0
  when /-j/
    job_name = arg.inspect.gsub('"', '')
  when /-v/
    var_column = arg.inspect.gsub('"', '').to_i
  when /-x/
    x_column = arg.inspect.gsub('"', '').to_i
  when /-l/
    take_log10 = true
  when /-o/
    output_name = arg.inspect.gsub('"', '')
  when /-s/
    slice_str = arg.inspect.gsub('"', '')
  end
end

puts "    job=#{job_name} x_column=#{x_column} var_column=#{var_column}"
puts "    output=#{output_name} slice_str=#{slice_str}"
if take_log10 then
  puts "    Will take logarithm of variable values."
end

tindx_file = "config/"+job_name+".times"
fin = File.open(tindx_file)
tindx_list = []; time_list = []
fin.each do |line|
  if line.start_with?("#") then next end
  if line.empty? then next end
  items = line.split()
  tindx_list << items[0].to_i
  time_list << items[1].to_f
end
fin.close
puts "Number of times found: #{tindx_list.count}"

fout = File.new(output_name, "w")
fout.puts "# x t var"

tindx_list.each do |tindx|
  cmd = "e4shared --post --job=#{job_name} --tindx-plot=#{tindx} --slice-list=#{slice_str}"
  o, e, s = Open3.capture3(*cmd.split)
  if s.success? then
    puts "tindx=#{tindx} time=#{time_list[tindx]} data-line-count=#{o.lines.count}"
    found_labels = false
    o.lines.each do |line|
      if line.match("pos.x") then
        found_labels = true
        labels = line.split
        # Remember that the labels line starts with a sharp character
        x_label = labels[x_column]
        var_label = labels[var_column]
        puts "    x_label=#{x_label} var_label=#{var_label}"
        next
      end
      if found_labels then
        if line.start_with?("Done") then break end
        items = line.split
        # The column labels start at 1
        x = items[x_column-1].to_f
        val = items[var_column-1].to_f
        if take_log10 then val = Math.log10(val) end
        fout.puts "#{x} #{time_list[tindx]} #{val}"
      end
    end
    fout.puts("")
  end
end

fout.close()

puts "Done."
