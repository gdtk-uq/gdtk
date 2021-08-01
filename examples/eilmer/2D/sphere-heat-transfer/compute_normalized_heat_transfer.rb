#! /usr/bin/env ruby
# compute_normalized_heat_transfer.rb
# Collect and plot the heat-transfer data from the loads files.
# PJ 2021-08-01

file = File.open("./loads/sphere-loads.times")
tindx = 0
file.each do |line|
  if line.start_with?('#') then next end
  items = line.split(' ')
  tindx = items[0].to_i
end
file.close
puts "tindx=#{tindx}"

values = {'x'=>Array.new, 'y'=>Array.new, 'q'=>Array.new}

(0..3).each do |i|
  filename = sprintf("./loads/t%04d/b%04d.t%04d.loads.dat", tindx, i, tindx)
  file = File.open(filename)
  columns = {'x'=>0, 'y'=>0, 'q'=>0}
  file.each do |line|
    if line.match?('pos.x')
      # pp line.match(/(\d+):pos.x\s/)
      columns['x'] = line.match(/(\d+):pos.x\s/)[1].to_i
      columns['y'] = line.match(/(\d+):pos.y\s/)[1].to_i
      columns['q'] = line.match(/(\d+):q_total\s/)[1].to_i
    end
    if line.start_with?('#') then next end
    if columns['x'] > 0
      items = line.strip.split(' ')
      values['x'] << items[columns['x']-1].to_f
      values['y'] << items[columns['y']-1].to_f
      values['q'] << items[columns['q']-1].to_f
    end
  end
  file.close
end

file = File.open("normalized_heat_transfer.data", 'w')
file.puts "# theta_degrees  q_W/m^^2  q/q_stag"
data = values['x'].zip(values['y'], values['q'])
data.each do |x, y, q|
  theta = Math.atan2(y,-x)*180.0/Math::PI
  q_norm = q/values['q'][0]
  file.puts "#{theta} #{q.abs} #{q_norm}"
end
file.close
