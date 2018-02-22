-- sketch-domain.lua
-- Called by the user input script to make a sketch of the flow domain.
-- PJ, 2016-08-13
s = Sketch:new{renderer="svg", projection="xyortho"}
s:set{canvas={0.0,0.0,120.0,80.0}, viewport={-2,-2,9,9}}

s:start{file_name="corner.svg"}
s:text{point=Vector3:new{x=4,y=8},
       text="Flow domain for backwards facing step",
       font_size=10}
s:text{point=Vector3:new{x=1,y=5},
       text="Blk0",
       font_size=10}
s:text{point=Vector3:new{x=4.5,y=5},
       text="Blk1",
       font_size=10}
s:text{point=Vector3:new{x=4.5,y=1.5},
       text="Blk2",
       font_size=10}

-- for flow domain
s:set{line_width=0.1, fill_colour="green"}
s:render{surf=quad0}; s:render{surf=quad1}; s:render{surf=quad2}

s:set{line_width=0.4}; s:render{path=ab}; s:render{path=cb} -- thick line for cone surface
s:dotlabel{point=a, label="a"}; s:dotlabel{point=b, label="b"}
s:dotlabel{point=c, label="c"}; s:dotlabel{point=d, label="d"}
s:dotlabel{point=e, label="e"}; s:dotlabel{point=f, label="f"}
s:dotlabel{point=g, label="g"}; s:dotlabel{point=h, label="h"}

s:set{line_width=0.3} -- for drawing rules
s:rule{direction="x", vmin=0.0, vmax=7.0, vtic=1,
       anchor_point=Vector3:new{x=0,y=-1},
       tic_mark_size=0.02, number_format="%.1f",
       text_offset=0.06, font_size=10}
s:text{point=Vector3:new{x=-0.5,y=-1}, text="x", font_size=12}
s:rule{direction="y", vmin=0.0, vmax=7.0, vtic=1,
       anchor_point=Vector3:new{x=-1,y=0},
       tic_mark_size=0.02, number_format="%.1f",
       text_offset=0.02, font_size=10}
s:text{point=Vector3:new{x=-1.0,y=-0.5}, text="y", font_size=12}

s:finish{}
