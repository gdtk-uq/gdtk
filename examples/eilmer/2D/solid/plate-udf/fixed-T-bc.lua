local T = {}
T[north] = 1000.0
T[east] = 800.0
T[west] = 300.0
T[south] = 300.0

function solidInterface(args)
   return T[whichBoundary]
end
