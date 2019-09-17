-- sketch-domain.lua
-- Called by the user input script to make a sketch of the flow domain.
-- PJ, 2018-03-17, 2018-09-03, 2019-09-03
--
-- During development of the geometry, I actually used renderer="xplot"
-- For making a document renderer="svg"
sk = Sketch:new{renderer="svg", projection="xyortho",
                canvas_mm={0.0,0.0,120.0,120.0}}
sk:set{viewport={-0.15,-0.15,0.2,0.2}}

sk:start{file_name="aerospike-domain.svg"}
sk:text{point=Vector3:new{x=0.025,y=0.18},
       text="flow domain for aerospike",
       font_size=16}

-- for flow domain, colour in patches
sk:set{line_width=0.1, fill_colour="green"}
for k,v in pairs(patches) do
   sk:render{surf=v}
end

-- thick line for nozzle surface and combustor walls
sk:set{line_width=0.5}
sk:render{path=lines.TD}
sk:render{path=lines.s0}
sk:render{path=lines.n0}
sk:render{path=lines.s3}
     
-- dotted labels for the defining points
some_points = {E=pnts.E, T=pnts.T, D=pnts.D, TDb1=pnts.TDb1, TDb2=pnts.TDb2}
for k,v in pairs(some_points) do
   sk:dotlabel{point=v, label=k}
end

-- a couple of axes, to give scale
sk:set{line_width=0.3} -- for drawing rules
sk:rule{direction="x", vmin=-0.05, vmax=0.15, vtic=0.05,
        anchor_point=Vector3:new{x=-0.05,y=-0.02},
        tic_mark_size=0.01, number_format="%.2f",
        text_offset=0.02, font_size=10}
sk:text{point=Vector3:new{x=0.05,y=-0.06}, text="x", font_size=12}
sk:rule{direction="y", vmin=0.0, vmax=0.15, vtic=0.05,
        anchor_point=Vector3:new{x=-0.080,y=0},
        tic_mark_size=0.01, number_format="%.2f",
        text_offset=0.015, font_size=10}
sk:text{point=Vector3:new{x=-0.125,y=0.075}, text="y", font_size=12}

sk:finish{}
