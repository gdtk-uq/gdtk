-- sketch-domain-extended.lua
-- Called by the user input script to make a sketch of the flow domain.
-- PJ, 2016-09-24
s = Sketch:new{renderer="svg", projection="xyortho"}
s:set{canvas={0.0,0.0,120.0,120.0}, viewport={-0.5,-0.5,2.5,2.5}}

s:start{file_name="conepe.svg"}
s:text{point=Vector3:new{x=1.0,y=2.2},
       text="extended flow domain",
       font_size=16}

-- for flow domain
s:set{line_width=0.1, fill_colour="green"}
s:render{surf=quad0}; s:render{surf=quad1}; s:render{surf=quad2}

s:set{line_width=0.5}; s:render{path=bc} -- thick line for cone surface
s:dotlabel{point=a, label="a"}; s:dotlabel{point=b, label="b"}
s:dotlabel{point=c, label="c"}; s:dotlabel{point=d, label="d"}
s:dotlabel{point=e, label="e"}; s:dotlabel{point=f, label="f"}

s:set{line_width=0.3} -- for drawing rules
s:rule{direction="x", vmin=0.0, vmax=2.0, vtic=0.5,
       anchor_point=Vector3:new{x=0,y=-0.1},
       tic_mark_size=0.05, number_format="%.1f",
       text_offset=0.12, font_size=10}
s:text{point=Vector3:new{x=1.0,y=-0.4}, text="x", font_size=12}
s:rule{direction="y", vmin=0.0, vmax=2.0, vtic=0.5,
       anchor_point=Vector3:new{x=-0.1,y=0},
       tic_mark_size=0.05, number_format="%.1f",
       text_offset=0.05, font_size=10}
s:text{point=Vector3:new{x=-0.35,y=1.0}, text="y", font_size=12}

s:finish{}
