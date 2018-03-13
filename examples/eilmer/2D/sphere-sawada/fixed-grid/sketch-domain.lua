-- sketch-domain.lua
-- Called by the user input script to make a sketch of the flow domain for Sawada's sphere.
-- PJ, 2016-08-31
s = Sketch:new{renderer="svg", projection="xyortho", canvas_mm={0.0,0.0,120.0,120.0}}
s:set{viewport={-0.090,-0.020,0.060,0.130}}

s:start{file_name="sphere.svg"}
s:text{point=Vector3:new{x=-0.010,y=0.110}, text="flow domain for sphere", font_size=16}

s:set{line_width=0.1, fill_colour="green"} -- for flow domain
s:render{surf=psurf}

s:set{line_width=0.5} -- thick line for sphere surface
s:render{path=bc}; s:render{path=cd}; s:render{path=de}
s:set{line_width=0.2} -- construction lines for Bezier curves
s:line{p0=i, p1=j, dashed=true}; s:line{p0=h, p1=k, dashed=true}
s:line{p0=e, p1=f, dashed=true}; s:line{p0=h, p1=g, dashed=true}

s:dotlabel{point=a, label="a"}; s:dotlabel{point=b, label="b"}
s:dotlabel{point=c, label="c"}; s:dotlabel{point=d, label="d"}
s:dotlabel{point=e, label="e"}; s:dotlabel{point=f, label="f"}
s:dotlabel{point=g, label="g"};
s:dotlabel{point=h, label="h"}; s:dotlabel{point=i, label="i"}
s:dotlabel{point=j, label="j"}; s:dotlabel{point=k, label="k"}

s:set{line_width=0.3} -- for drawing rules
s:rule{direction="x", vmin=-0.060, vmax=0.040, vtic=0.020, anchor_point=Vector3:new{x=0,y=-0.002},
       tic_mark_size=0.001, number_format="%.3f", text_offset=0.007, font_size=10}
s:text{point=Vector3:new{x=-0.010,y=-0.015}, text="x,m", font_size=12}
s:rule{direction="y", vmin=0.0, vmax=0.120, vtic=0.020, anchor_point=Vector3:new{x=-0.060,y=0},
       tic_mark_size=0.001, number_format="%.3f", text_offset=0.002, font_size=10}
s:text{point=Vector3:new{x=-0.075,y=0.070}, text="y,m", font_size=12}

s:finish{}
