#import "lib.typ": cram-snap, theader
//#import "@preview/cram-snap:0.2.2": cram-snap, theader


#set page(
  paper: "a4",
  flipped: true,
  margin: 1cm,
//  background: rotate(-45deg,text(font:"Arial",size:200pt, fill: rgb("FFEEEE"))[DRAFT])
)
//#set text(font: "Courier", size: 11pt)

#show: cram-snap.with(
  title: [EILMER CHEATSHEET _(v5.0.0)_],
  subtitle: [Common commands and files],
  icon: image("eilmer-icon-2.png"),
)

#set table(columns: (1.2fr, 1fr), align: (top))

#table(
  theader[Meta commands ],
  [`lmr help`], [print general help],
  [`lmr help -a`], [list all commands],
  [`lmr help <cmd-name>`], [print help about _cmd-name_],
  [`lmr version`], [print version information related to compilation],
  [`lmr revision-id`], [print repository revision ID]
)  

#table(
  theader[Preparation],
  [`lmr prep-gas -i <input> -o <output>`], [prepare a gas model file],
  [`lmr prep-reactions -g <gasfile>` \
   `    -i <input> -o <output>`],[prepare reaction scheme],
  [`lmr prep-grid [-j <file>]`], [prepare a grid (default: `job.lua`)],
  [`lmr prep-sim  [-j <file>]`], [prepare a simulation (default: `job.lua`)],
)

#table(
  theader[Running simulations],
  [`lmr run`], [run using shared memory],
  [`lmr run -s 2`], [restart simulation from snapshot 2],
  [`lmr run --max-wall-clock=<s>`], [limit run time to _s_ seconds],
  [`mpirun -np 8 lmr-mpi-run`], [run MPI code using real-valued numbers],
  [`mpirun -np 8 lmrZ-mpi-run`], [run MPI code using complex-valued numbers],
  [`lmr plot-diagnostics`], [real-time or offline residual monitoring]
)

#colbreak()

#table(
  theader[Post processing],
  [`lmr snapshot2vtk`], [convert field data to VTK format],
  [`lmr snapshot2vtk --add-vars="mach"`], [add field variables when generating VTK],
  [`lmr slice-flow -n "var1,var2" -l "blk-range,i-range,j-range,k-range" -o slice.dat`], [extract variables from a slice through the domain, place in `slice.dat` file],
  [`lmr extract-line -l "x0,y0,z0,x1,y1,z1,n" -o line.dat`], [extract line of data using _n_ samples from $(x_0, y_0, z_0)$ to $(x_1, y_1, z_1)$, place in `line.dat` file],
  [`lmr probe-flow -l "x0,y0,z0"` \ `[-o probe.dat]`], [sample flow field at location extract $(x, y, z)$ [, place in `probe.dat` file]],
  
)

#table(fill: none, 
  theader[lmr files overview],
  table.cell(colspan: 2, align: center)[#set text(font: "Roboto Mono"); #image("simulation-process-chart.svg", width: 85%)],
)


