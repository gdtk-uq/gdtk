-- sketch-domain.lua
-- Called by the user input script to make a sketch of the flow domain for Sawada's sphere.
-- PJ, 2016-08-31
s = Sketch:new{renderer="svg", projection="xyortho", canvas_mm={0.0,0.0,120.0,120.0}}
s:set{viewport={-0.090,-0.020,0.060,0.130}}

s:start{file_name="sphere.svg"}
s:text{point=Vector3:new{x=-0.010,y=0.110}, text="flow domain for sphere", font_size=16}

s:set{line_width=0.1, fill_colour="green"} -- for flow domain
s:render{surf=bp.patch}

s:set{line_width=0.3} -- for drawing rules
s:rule{direction="x", vmin=-0.060, vmax=0.040, vtic=0.020, anchor_point=Vector3:new{x=0,y=-0.002},
       tic_mark_size=0.001, number_format="%.3f", text_offset=0.007, font_size=10}
s:text{point=Vector3:new{x=-0.010,y=-0.015}, text="x,m", font_size=12}
s:rule{direction="y", vmin=0.0, vmax=0.120, vtic=0.020, anchor_point=Vector3:new{x=-0.060,y=0},
       tic_mark_size=0.001, number_format="%.3f", text_offset=0.002, font_size=10}
s:text{point=Vector3:new{x=-0.075,y=0.070}, text="y,m", font_size=12}

s:finish{}
