#!/usr/bin/env ruby
# Prepare to restart from an eilmer snapshot by placing
# in an appropriate location in the time series.
#
# Author: Rowan J. Gollan
# Version: 2020-04-23 -- first go

require 'getoptlong'
require 'open3'
require 'fileutils'

opts = GetoptLong.new(
  ["--job", "-j", GetoptLong::REQUIRED_ARGUMENT],
  ["--snapshot", "-s", GetoptLong::REQUIRED_ARGUMENT],
  ["--verbose", "-v", GetoptLong::NO_ARGUMENT],
  ["--replace", "-r", GetoptLong::OPTIONAL_ARGUMENT],
  ["--help", "-h", GetoptLong::NO_ARGUMENT]
)

puts "e4restart"
puts "Prepare an eilmer restart from a snapshot."

jobName = ""
snapshotIdx = -1
replaceIdx = -1
$verbose = false

opts.each do |opt, arg|
  case opt
  when /-h/
    puts "Usage:"
    puts "$ e4-prep-restart --job=<jobName> --snapshot=<s_index> ?--replace=<t_index>? ?--verbose?"
    puts "Notes:"
    puts "We work with the snapshot specified by <s_index>."
    puts "The default action of this program is append the snapshot to the time series."
    puts "If, instead, one wants to replace a solution in the time series with the snapshot,"
    puts "the replace option (--replace) is available."
    puts "Option --replace=<t_index>|-r <t_index> replaces the a solution in the time series at <t_index> with snapshot."
    puts "Option --verbose|-v be a bit noisy about what files are being moved and copied."
    exit 0
  when /-j/
    jobName = arg.inspect.gsub('"', '')
  when /-s/
    snapshotIdx = arg.to_i
  when /-r/
    replaceIdx = arg.to_i
    if replaceIdx < 0
      puts "ERROR: Invalid replace index provided: --replace=#{i}"
      exit 1
    end
  when /-v/
    $verbose = true
  end
end

if snapshotIdx == -1
  puts "ERROR: No snapshot index supplied. Use --snapshot=<s_index>"
  exit 1
end

if jobName == ""
  puts "ERROR: No <jobName> suplied. Use --job=my_job"
  exit 1
end

def readTimesFile(fname)
  # We store the read in times index as a hash in case there are duplicate entries
  tInfo = Hash.new
  lines = IO.readlines(fname)
  lines.drop(1).each do |line|
    tks = line.split(" ")
    tInfo[tks[0].to_i] = {"t" => tks[1].to_f, "dt" => tks[2].to_f}
  end
  return tInfo.sort, tInfo
end

def writeTimesFile(fname, tInfo)
  File.open(fname, "w") do |f|
    f.puts "# tindx sim_time dt_global"
    tInfo.each_with_index do |inf, idx|
      f.puts "%04d %.18e %.18e" % [idx, inf["t"], inf["dt"]]
    end
  end
end

def copyFiles(jobName, type, snapshotIdx, tIdx)
  FileUtils.mkdir_p("#{type}/t%04d" % [tIdx], verbose: $verbose)
  files = Dir.glob("#{type}/snapshot-%04d/#{jobName}*" % [snapshotIdx])
  files.each do |f|
    tks = f.split("/")
    nps = tks[2].split(".")
    FileUtils.cp(f, "#{type}/t%04d/#{nps[0]}.#{nps[1]}.#{nps[2]}.t%04d.#{nps[5]}" % [tIdx, tIdx], verbose: $verbose)
  end
end

tgtTypes = ["flow", "solid"]

if replaceIdx >= 0
  # We insert the snapshot into the time series
  # The two actions are to:
  #   1. Copy the appropriate snapshot files across
  #   2. Insert the new time,dt info into the .times file
  tFile = "config/#{jobName}.times"
  timesInfo, tHash = readTimesFile(tFile)
  # 1. Copy snapshot into replacement position
  # check that the entry we want to replace does exist
  if !tHash.key?(replaceIdx)
    puts "ERROR: The time series index to be replaced could not be found: #{replaceIdx}"
    exit 1
  end
  snapInfo,_ = readTimesFile("config/#{jobName}.snapshots")
  if snapshotIdx >= snapInfo.count
    puts "ERROR: The supplied snapshot index #{snapshotIdx} is not available in the series."
    exit 1
  end
  # Take a peek to see if there is a time series of grids.
  # If so, we add that to the targets for copying across.
  tgtTypes << "grid" if Dir.exist?("grid/t%04d" % [replaceIdx])
  # Gather up flow files and copy over
  tgtTypes.each do |type| 
    copyFiles(jobName, type, snapshotIdx, replaceIdx)
  end
  # 2. Alter .times file to include information from the replacment snapshot
  # We're going to work in the file from reverse because sometimes
  # there are multiple duplicate time index entries if there has not
  # been clean starts.
  lines = IO.readlines(tFile)
  lineToReplace = 0
  lines.to_enum.with_index.reverse_each do |line, i|
    tks = line.split(" ")
    if tks[0].to_i == replaceIdx
      lineToReplace = i + 1 # because sed counts lines from 1
      break
    end
  end
  # We'll use sed to replace the line since this seems the simplest
  snap = snapInfo[snapshotIdx][1]
  sedCmd = "sed -i '#{lineToReplace} s/.*/%04d %.18e %.18e/' #{tFile}" % [replaceIdx, snap["t"], snap["dt"]]
  print sedCmd if $verbose
  o, e, s = Open3.capture3(sedCmd)
  if s.success?
    puts " : ok" if $verbose
  else
    puts " : fail" if $verbose
    raise "sed error: #{e}"
  end
  puts "Snapshot has been added at tindx: #{replaceIdx}"
else
  # We append the snapshot to the time series
  # The two actions are to:
  #   1. Copy the appropriate snapshot files across
  #   2. Append appropriate information to the times file
  timesInfo,_ = readTimesFile("config/#{jobName}.times")
  snapInfo,_ = readTimesFile("config/#{jobName}.snapshots")
  if snapshotIdx >= snapInfo.count
    puts "ERROR: The supplied snapshot index #{snapshotIdx} is not available in the series."
    exit 1
  end
  # 1. Copy snapshot into next time series slot
  tidxNew = timesInfo[-1][0] + 1
  # Take a peek to see if there is a time series of grids.
  # If so, we add that to the targets for copying across.
  tgtTypes << "grid" if Dir.exist?("grid/t%04d" % [tidxNew-1])
  # Gather up flow files and copy over
  tgtTypes.each do |type| 
    copyFiles(jobName, type, snapshotIdx, tidxNew)
  end
  # 2. Append to the times file
  snap = snapInfo[snapshotIdx][1]
  File.open("config/#{jobName}.times", "a") do |f|
    f.puts "%04d %.18e %.18e" % [tidxNew, snap["t"], snap["dt"]]
  end
  puts "Snapshot has been added at tindx: #{tidxNew}"
end



                     
