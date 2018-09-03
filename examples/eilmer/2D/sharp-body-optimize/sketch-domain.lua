-- sketch-domain.lua
-- Called by the user input script to make a sketch of the flow domain.
-- PJ, 2018-03-17, 2018-09-03
--
-- During development of the geometry, I actually used renderer="xplot"
-- For making a document renderer="svg"
sk = Sketch:new{renderer="svg", projection="xyortho",
                canvas_mm={0.0,0.0,120.0,120.0}}
sk:set{viewport={-2.5,-1.0,7.0,8.5}}

sk:start{file_name="sharp-body.svg"}
sk:text{point=Vector3:new{x=2.5,y=8.0},
       text="flow domain for sharp body",
       font_size=16}

-- for flow domain, colour in patches
sk:set{line_width=0.1, fill_colour="green"}
for k,v in pairs(patches) do
   sk:render{surf=v}
end

-- thick line for plate surface
sk:set{line_width=0.5}
sk:render{path=lines.body}
     
-- dotted labels for all of the defining points
for k,v in pairs(pnts) do
   sk:dotlabel{point=v, label=k}
end

-- a couple of axes, to give scale
sk:set{line_width=0.3} -- for drawing rules
sk:rule{direction="x", vmin=-1.0, vmax=6.0, vtic=1.0,
        anchor_point=Vector3:new{x=-1.0,y=-0.2},
        tic_mark_size=0.2, number_format="%.2f",
        text_offset=0.5, font_size=10}
sk:text{point=Vector3:new{x=2.5,y=-1.0}, text="x", font_size=12}
sk:rule{direction="y", vmin=0.0, vmax=8.0, vtic=1.0,
        anchor_point=Vector3:new{x=-1.2,y=0},
        tic_mark_size=0.2, number_format="%.2f",
        text_offset=0.3, font_size=10}
sk:text{point=Vector3:new{x=-2.3,y=4.5}, text="y", font_size=12}

sk:finish{}
