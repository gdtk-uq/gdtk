#!/usr/bin/env ruby
# Compact the many solution files for a job into a few archive files.
#
# Peter J. UQ
# Version: 2019-06-05 first experiment
#
require 'getoptlong'
require 'open3'
require 'fileutils'

opts = GetoptLong.new(
  ["--job", "-j", GetoptLong::REQUIRED_ARGUMENT],
  ["--verbose", "-v", GetoptLong::NO_ARGUMENT],
  ["--force", "-f", GetoptLong::NO_ARGUMENT],
  ["--restore", "-r", GetoptLong::NO_ARGUMENT],
  ["--help", "-h", GetoptLong::NO_ARGUMENT]
)

puts "e4compact"
puts "Compact the many solution files for a job into a few archive files."

job_name = ""
target_dirs = ["grid", "flow", "solid", "loads", "hist", "plot", "CellData"]
$verbose = false
$force = false
$restore = false

opts.each do |opt, arg|
  case opt
  when /-h/
    puts "Usage:"
    puts "$ e4compact --job=<jobName> ?--verbose? ?--force? ?--restore?"
    puts "Notes:"
    puts "We delegate most of the work to the tar command."
    puts "Archive files will have names of the form jobName-dirName.tar."
    puts "Solution files that are successfully archived are removed from the directories."
    puts "If the directories are then empty, they are also removed."
    puts "Existing archives will only be updated if --force|-f option is specified."
    puts "Option --restore|-r will extract files from archives and restore the directories."
    puts "Option --verbose|-v will echo the tar commands as they are issued."
    exit 0
  when /-j/
    job_name = arg.inspect.gsub('"', '')
  when /-v/
    $verbose = true
  when /-f/
    $force = true
  when /-r/
    $restore = true
  end
end

def compact_dir(dir_name, job_name)
  #
  # First, consider if we should compact the directory.
  if !Dir.exist?(dir_name)
    puts "Do not see directory #{dir_name}"
    return
  end
  archive_file = job_name+"-"+dir_name+".tar"
  puts "Compact directory #{dir_name} into archive #{archive_file}"
  if File.exist?(archive_file) && !$force
    puts "  Archive already exists. Will not update."
    return
  end
  #
  # Gather up the files in the directory; there may be many.
  case dir_name
  when "grid"
    files = Dir.glob("#{dir_name}/t????/#{job_name}*")
  when "flow"
    files = Dir.glob("#{dir_name}/t????/#{job_name}*")
  when "solid"
    files = Dir.glob("#{dir_name}/t????/#{job_name}*")
  when "loads"
    files = Dir.glob("#{dir_name}/t????/b*.dat")
    files += Dir.glob("#{dir_name}/#{job_name}-loads.times")
  when "hist"
    files = Dir.glob("#{dir_name}/#{job_name}*")
  when "plot"
    files = Dir.glob("#{dir_name}/#{job_name}*")
  else
    files = []
  end
  if files.length == 0
    puts "  No files found for job=#{job_name}. Delete directory if empty."
    FileUtils.rm_r(dir_name)
    return
  end
  #
  # Build the archive file at this level.
  files.each_with_index do |f, i|
    cmd = "tar --verify --file=#{archive_file}"
    if i == 0 && !File.exist?(archive_file)
      cmd << " --create"
    else
      cmd << " --update"
    end
    cmd << " #{f}"
    print "  cmd= \"#{cmd}\"" if $verbose
    o, e, s = Open3.capture3(*cmd.split)
    if s.success?
      puts " ok" if $verbose
    else
      puts " fail" if $verbose
      raise "tar error: #{e}"
    end
  end
  #
  # Since the archive was successfully built,
  # remove the files and directory, if it is empty.
  files.each do |f| File.delete(f) end
  FileUtils.rm_r(dir_name)
end

def restore_dir(dir_name, job_name)
  # Extract files from the archive.
  archive_file = job_name+"-"+dir_name+".tar"
  puts "Restore directory #{dir_name} from archive #{archive_file}"
  if File.exist?(archive_file)
    cmd = "tar -xf #{archive_file}"
    print "  cmd= \"#{cmd}\"" if $verbose
    o, e, s = Open3.capture3(*cmd.split)
    if s.success?
      puts " ok" if $verbose
    else
      puts " fail" if $verbose
      raise "tar error: #{e}"
    end
  end
end

input_script = job_name+".lua"
if File.exist?(input_script)
  puts "Can see job input script so do some work..."
  target_dirs.each do |name|
    if $restore
      restore_dir(name, job_name)
    else
      compact_dir(name, job_name)
    end
  end
else
  puts "Cannot see job input script: #{input_script}"
end

puts "Done."
