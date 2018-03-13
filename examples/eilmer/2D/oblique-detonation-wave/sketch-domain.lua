-- sketch-domain.lua
-- Called by the user input script to make a sketch of the flow domain.
-- PJ, 2017-01-07, adapted from cone20 case
s = Sketch:new{renderer="svg", projection="xyortho", canvas_mm={0.0,0.0,120.0,120.0}}
s:set{viewport={-0.7,-0.4,2.2,2.5}}

s:start{file_name="odw.svg"}
s:text{point=Vector3:new{x=0.75,y=2.2},
       text="oblique detonation wave domain",
       font_size=16}

-- for flow domain
s:set{line_width=0.1, fill_colour="green"}
s:render{surf=patch0}; s:render{surf=patch1}

s:set{line_width=0.5}; s:render{path=south1} -- thick line for wedge surface
s:dotlabel{point=a, label="a"}; s:dotlabel{point=b, label="b"}
s:dotlabel{point=c, label="c"}; s:dotlabel{point=d, label="d"}
s:dotlabel{point=e, label="e"}; s:dotlabel{point=f, label="f"}

s:set{line_width=0.3} -- for drawing rules
s:rule{direction="x", vmin=-0.5, vmax=2.0, vtic=0.5,
       anchor_point=Vector3:new{x=0,y=-0.1},
       tic_mark_size=0.02, number_format="%.1f",
       text_offset=0.12, font_size=10}
s:text{point=Vector3:new{x=0.75,y=-0.35}, text="x", font_size=12}
s:rule{direction="y", vmin=0.0, vmax=2.0, vtic=0.5,
       anchor_point=Vector3:new{x=-0.35,y=0},
       tic_mark_size=0.02, number_format="%.1f",
       text_offset=0.02, font_size=10}
s:text{point=Vector3:new{x=-0.5,y=1.25}, text="y", font_size=12}

s:finish{}
