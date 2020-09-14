-- sketch-domain.lua
-- Called by the user input script to make a sketch of the flow domain.
-- PJ, 2018-03-17, 2018-09-03, 2019-09-03, 2020-08-01
--
-- During development of the geometry, I actually used renderer="xplot"
-- For making a document renderer="svg"
sk = Sketch:new{renderer="svg", projection="xyortho",
                canvas_mm={0.0,0.0,120.0,120.0}}
sk:set{viewport={-0.20,-0.20,1.60,1.60}}

sk:start{file_name="scramjet-domain.svg"}
sk:text{point=Vector3:new{x=0.7,y=1.2},
       text="flow domain for scramjet inlet",
       font_size=16}

-- for flow domain, colour in patches
sk:set{line_width=0.1, fill_colour="green"}
for k,v in pairs(patches) do
   sk:render{surf=v}
end

-- thick line for compression surface and cowl
sk:set{line_width=0.5}
sk:render{path=lines.s0}
sk:render{path=lines.s1}
sk:render{path=lines.s2}
sk:render{path=lines.n1}
sk:render{path=lines.n2}

-- dotted labels for the defining points
for k,v in pairs(pnts) do
   sk:dotlabel{point=v, label=k}
end

-- a couple of axes, to give scale
sk:set{line_width=0.3} -- for drawing rules
sk:rule{direction="x", vmin=0.0, vmax=1.5, vtic=0.5,
        anchor_point=Vector3:new{x=-0.05,y=-0.02},
        tic_mark_size=0.03, number_format="%.2f",
        text_offset=0.08, font_size=10}
sk:text{point=Vector3:new{x=0.8,y=-0.15}, text="x", font_size=12}
sk:rule{direction="y", vmin=0.0, vmax=1.5, vtic=0.5,
        anchor_point=Vector3:new{x=-0.05,y=0},
        tic_mark_size=0.03, number_format="%.2f",
        text_offset=0.03, font_size=10}
sk:text{point=Vector3:new{x=-0.15,y=0.80}, text="y", font_size=12}

sk:finish{}
