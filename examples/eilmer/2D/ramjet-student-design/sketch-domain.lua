-- sketch-domain.lua
-- Called by the user input script to make a sketch of the flow domain.
-- PJ, 2016-09-26
sk = Sketch:new{renderer="svg", projection="xyortho", canvas_mm={0.0,0.0,120.0,120.0}}
sk:set{viewport={-0.2,-0.3,0.7,0.6}}

sk:start{file_name="ram2.svg"}
sk:text{point=Vector3:new{x=0.25,y=0.35},
	text="flow domain for student ramjet",
	font_size=16}

-- for flow domain
sk:set{line_width=0.1, fill_colour="green"}
for ib = 0, 15 do
   sk:render{surf=quad[ib]}
end

sk:set{line_width=0.5}; -- sk:render{path=bc} -- thick line for cone surface

-- sk:dotlabel{point=a, label="a"}; sk:dotlabel{point=b, label="b"}
sk:dotlabel{point=c, label="c"}; -- sk:dotlabel{point=d, label="d"}

sk:set{line_width=0.3} -- for drawing rules
sk:rule{direction="x", vmin=-0.1, vmax=0.6, vtic=0.1,
	anchor_point=Vector3:new{x=0,y=-0.21},
	tic_mark_size=0.01, number_format="%.1f",
	text_offset=0.03, font_size=10}
sk:text{point=Vector3:new{x=0.25,y=-0.27}, text="x", font_size=12}
sk:rule{direction="y", vmin=-0.2, vmax=0.30, vtic=0.1,
	anchor_point=Vector3:new{x=-0.10,y=0},
	tic_mark_size=0.01, number_format="%.1f",
	text_offset=0.02, font_size=10}
sk:text{point=Vector3:new{x=-0.17,y=0.05}, text="y", font_size=12}

sk:finish{}
