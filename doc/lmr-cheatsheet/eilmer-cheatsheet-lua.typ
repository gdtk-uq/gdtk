
#import "lib.typ": cram-snap, theader

#set page(
  paper: "a4",
  flipped: true,
  margin: 1cm,
  //background: rotate(-45deg,text(font:"Arial",size:200pt, fill: rgb("FFEEEE"))[DRAFT])
)

#show: cram-snap.with(
  title: [EILMER CHEATSHEET _(v5.0.0)_],
  subtitle: [Common Lua input],
  icon: image("eilmer-icon-2.png"),
)


#set table(columns: (0.6fr, 1fr), align: (top))

#table(
  theader[Common config settings],
  table.cell(colspan: 2)[#text(style: "oblique")[All options are accessed as `config.`_option_, e.g. `config.solver_mode`]],
  [`solver_mode`], [`transient` | `steady`],
  [`viscous`], [set `true` to include viscous terms (Navier-Stokes)],
  [`dimensions`], [`2` | `3`],
  [`axisymmetric`], [set `true` for 2D geometries with axial symmetry about $y = 0$],
  [`gasdynamic_update_scheme`], [common explicit schemes: \ `'predictor-corrector'`, `'classic-rk3'` \
                                        moving grid choices:  `'moving-grid-1-stage'`, `'moving-grid-2-stage'` \
                                        implicit schemes: \ `'backward-euler'`, `'implicit_rk1'`],
  [`cfl_value`], [ _float_ value for target Courant-Friedrichs-Lewy number],
  [`cfl_schedule`], [table of `{time, CFL}` pairs such as \ `{{0.0, 0.25}, {1.0e-4, 0.75}}`],
  [`dt_init`], [ _float_ for initial timestep],
  [`max_time`], [ _float_ for maximum simulated time],
  [`max_step`], [ _int_ for maximum number of update steps],
  [`dt_plot`], [_float_ period between writing out flow field],
  [`dt_history`], [_float_ period between writing out history data],
  [`interpolation_order`], [`2`: high order reconstruction \
                            `1`: no reconstruction],
  [`flux_calculator`], [common choices include: `'ausmdv'`, `'adaptive_hlle_asumdv'`, `'adaptive_hanel_ausmdv'`],
  
  
)

#colbreak()

#set table(columns: (1.2fr, 1fr))



#table(
  theader[Geometry and patches and clustering (Oh my!)],
  [`Vector3:new{x=x0, y=y0, z=z0}`], [Creates Cartesian vector object at $(x_0, y_0, z_0)$],
  [`Line:new{p0=a, p1=b}`], [ #line(length: 4cm, angle: +4deg)
                              #place(dx: 0pt, dy: -9.5pt, circle(fill: red, radius: 2pt))
                              #place(dx: 3.9cm, dy: -2pt, circle(fill: red, radius: 2pt))
                              #place(dx: -8pt, dy: -10.5pt, text[$arrow(a)$]) 
                              #place(dx: 4.1cm, dy: -2.5pt, text[$arrow(b)$]) ],
  [`Arc:new{p0=a, p1=b, centre=c}`], [Arc from $arrow(a)$ to $arrow(b)$ about centre $arrow(c)$],
  [`Arc3:new{p0=a, pmid=b, p1=c}` \ ` ` ], [Arc through 3 points: $arrow(a)$, $arrow(b)$ and $arrow(c)$ \ ` `],
  [#text(size: 10pt)[`Bezier:new{points={p0, p1, p2, p3}}`]], [ #place(dx: 0.0cm - 1pt, dy: 0.5cm - 1pt, circle(fill: orange, radius: 2pt))
                                             #place(dx: 1.5cm - 1pt, dy: -0.3cm - 1pt, circle(fill: orange, radius: 2pt))
                                             #place(dx: 3cm - 1pt, dy: 0.7cm - 1pt, circle(fill: orange, radius: 2pt))
                                             #place(dx: 5cm - 1pt, dy: 0cm - 2pt, circle(fill: orange, radius: 2pt))
                                             #curve(stroke: black, curve.move((0cm,+0.5cm)), curve.cubic((1.5cm, -0.3cm), (3cm, +0.7cm), (5cm, 0cm)))
                                             #place(dx: 0.0cm - 3pt, dy: +5pt, text[$arrow(p)_0$])
                                             #place(dx: 1.7cm, dy: -25pt, text[$arrow(p)_1$])
                                             #place(dx: 3.2cm, dy: +8pt, text[$arrow(p)_2$])
                                             #place(dx: 5.2cm, dy: -15pt, text[$arrow(p)_3$])
                                            ],
  [ \ `Spline:new{points={p0, p1, ...}}`], [ \ Spline through points $(p_0, p_1, ...)$],
  [ `Spline2:new{filename='file.txt'}`], [ Spline through points in `file.txt`],
  [ `Polyline:new{segments={pathA, pathB, ...}}`], [ Single path built from joined segments],
  [ `CoonsPatch:new{north=..., south=..., east=..., west=...}`], [ Coons patch with edges _N_, _S_, _E_ and _W_],
  [ `AOPatch:new{north=..., south=..., east=..., west=..., [nx=..., ny=...]}`], [Area-orthogonal patch with background grid dimensions $(n_x, n_y)$],
  [ `RobertsFunction:new{end0, end1, beta}`], [Roberts clustering function. Set `end0|1` true to cluster towards `end0|1`.
  `beta` controls strength of clustering. ],
  [ `GaussianFunction:new{m, s, ratio}`], [Gaussian clustering centred at $m$ with width $s$ and tightness controlled by _ratio_],
    
)

#table(
  theader[Grids and block construction],
  [`StructuredGrid:new{psurface, niv, njv, cfList}`], [Initialise `StructuredGrid` object],
  [`UnstructuredGrid:new{sgrid}`], [Initialise `UnstructuredGrid` from a structured grid],
  [ #{show raw: set text(size: 8pt); `UnstructuredGrid:new{filename="grid.su2", fmt="su2text"}`}], [Initialise `UnstructuredGrid` from `grid.su2` in SU#super[2] format],
  [`registerFluidGrid{grid=..., fsTag=..., bcTags=...}`], [Register grids at `prep-grid` stage],
  [`identifyGridConnections()`], [Automate search for grid connections and apply],
  [`makeFluidBlocks(bcDict, flowDict)`], [From registered grids, make fluid blocks],
   
)
  



