-- sketch-domain.lua
-- Called by the user input script to make a sketch of the flow domain for sphere-heat-transfer.
-- PJ, 2016-08-31
s = Sketch:new{renderer="svg", projection="xyortho"}
s:set{canvas={0.0,0.0,120.0,120.0}, viewport={-0.020,-0.005,0.010,0.025}}

s:start{file_name="sphere.svg"}
s:text{point=Vector3:new{x=-0.007,y=0.022}, text="flow domain for sphere", font_size=16}

s:set{line_width=0.1, fill_colour="green"} -- for flow domain
s:render{surf=psurf}

s:set{line_width=0.5}; s:render{path=sphere_edge} -- thick line for sphere surface
s:dotlabel{point=a, label="a"}; s:dotlabel{point=b, label="b"}
s:dotlabel{point=c, label="c"}
-- points defining the Bezier curve along the inflow boundary
for i=1, #d do
   s:dotlabel{point=d[i], label=string.format("d[%d]", i)}
end

s:set{line_width=0.3} -- for drawing rules
s:rule{direction="x", vmin=-0.015, vmax=0.0, vtic=0.005, anchor_point=Vector3:new{x=0,y=-0.001},
       tic_mark_size=0.0002, number_format="%.3f", text_offset=0.0012, font_size=10}
s:text{point=Vector3:new{x=-0.007,y=-0.0035}, text="x,m", font_size=12}
s:rule{direction="y", vmin=0.0, vmax=0.020, vtic=0.005, anchor_point=Vector3:new{x=-0.013,y=0},
       tic_mark_size=0.0002, number_format="%.3f", text_offset=0.001, font_size=10}
s:text{point=Vector3:new{x=-0.018,y=0.0125}, text="y,m", font_size=12}

s:finish{}
