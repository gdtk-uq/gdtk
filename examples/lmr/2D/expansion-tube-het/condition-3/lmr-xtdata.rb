#!/usr/bin/env ruby
# lmr-xtdata.rb
# Accumulate x,t data for the selected variable, over the tindx saved solutions.
# PJ, 2019-10-22 initial experiment for Hans H. and David G.
#     2019-10-23 Make the script a bit more versatile.
#     2025-02-10 Port to lmr(5).
#
require 'getoptlong'
require 'open3'
require 'fileutils'
require 'yaml'

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

# Default values for parameters
x_column = "pos.x"
var_column = "p"
take_log10 = false
output_name = "xt.data"
slice_str = ":,:,0,0"

opts.each do |opt, arg|
  case opt
  when /-h/
    puts "Usage:"
    puts "$ xtdata.rb ?--vcolumn=p? ?--xcolumn=pos.x?"
    puts "    ?--log10? ?--output=xt.data? ?--slice=:,:,0,0?"
    puts "The default slice string assumes that all blocks are in a single line"
    puts "along the x-axis and, for each block, you want to select the full i-range"
    puts "of cells with j=0 and k=0 cell indices."
    exit 0
  when /-j/
    job_name = arg.inspect.gsub('"', '')
  when /-v/
    var_column = arg.inspect.gsub('"', '')
  when /-x/
    x_column = arg.inspect.gsub('"', '')
  when /-l/
    take_log10 = true
  when /-o/
    output_name = arg.inspect.gsub('"', '')
  when /-s/
    slice_str = arg.inspect.gsub('"', '')
  end
end

puts "    x_column=#{x_column} var_column=#{var_column}"
puts "    output=#{output_name} slice_str=#{slice_str}"
if take_log10 then
  puts "    Will take logarithm of variable values."
end

times_file = "./lmrsim/snapshots/snapshot-times-metadata"
times_metadata = YAML.load_file(times_file)
snaps = times_metadata.keys.sort
puts "Number of snapshots found: #{snaps.count}"

fout = File.new(output_name, "w")
fout.puts "# x t var"

snaps.each do |snap|
  cmd = "lmr slice-flow --snapshot=#{snap} --slice-list=#{slice_str} --names=#{var_column}"
  puts "cmd= #{cmd}"
  o, e, s = Open3.capture3(*cmd.split)
  if s.success? then
    t = times_metadata[snap]["time"]
    puts "snap=#{snap} time=#{t} data-line-count=#{o.lines.count}"
    found_labels = false
    x_index = -1
    var_index = -1
    o.lines.each do |line|
      if line.match("pos.x") then
        found_labels = true
        labels = line.split
        x_index= labels.index(x_column)
        var_index= labels.index(var_column)
        puts "    x_index=#{x_index} var_index=#{var_index}"
        next
      end
      if found_labels then
        if line.start_with?("Done") then break end
        items = line.split
        # The column labels start at 1
        x = items[x_index].to_f
        val = items[var_index].to_f
        if take_log10 then val = Math.log10(val) end
        fout.puts "#{x} #{t} #{val}"
      end
    end
    fout.puts("")
  end
end

fout.close()

puts "Done."
