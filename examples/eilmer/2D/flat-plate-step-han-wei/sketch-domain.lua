-- sketch-domain.lua
-- Called by the user input script to make a sketch of the flow domain.
-- PJ, 2018-03-17
--
-- During development of the geometry, I actually used renderer="xplot"
sk = Sketch:new{renderer="svg", projection="xyortho",
                canvas_mm={0.0,0.0,120.0,120.0}}
sk:set{viewport={-0.03,-0.03,0.15,0.15}}

sk:start{file_name="plate.svg"}
sk:text{point=Vector3:new{x=0.06,y=0.06},
       text="Han Wei's heated plate",
       font_size=16}

-- for flow domain, colour in quadrilaterals
sk:set{line_width=0.1, fill_colour="green"}
for k,v in pairs(quads) do
   sk:render{surf=v}
end

-- thick line for plate surface
sk:set{line_width=0.5}
sk:render{path=Line:new{p0=pnts.a,p1=pnts.e}}
sk:render{path=Line:new{p0=pnts.e,p1=pnts.j}}
sk:render{path=Line:new{p0=pnts.j,p1=pnts.k}}
     
-- dotted labels for all of the defining points
for k,v in pairs(pnts) do
   sk:dotlabel{point=v, label=k}
end
-- a couple of axes, to give scale
sk:set{line_width=0.3} -- for drawing rules
sk:rule{direction="x", vmin=0.0, vmax=0.12, vtic=0.02,
        anchor_point=Vector3:new{x=0,y=-0.005},
        tic_mark_size=0.002, number_format="%.2f",
        text_offset=0.006, font_size=10}
sk:text{point=Vector3:new{x=0.07,y=-0.015}, text="x", font_size=12}
sk:rule{direction="y", vmin=0.0, vmax=0.04, vtic=0.02,
        anchor_point=Vector3:new{x=-0.005,y=0},
        tic_mark_size=0.002, number_format="%.2f",
        text_offset=0.002, font_size=10}
sk:text{point=Vector3:new{x=-0.02,y=0.03}, text="y", font_size=12}

sk:finish{}
