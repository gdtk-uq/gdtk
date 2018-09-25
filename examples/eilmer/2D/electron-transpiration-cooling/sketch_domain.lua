-- sketch-domain.lua
-- Called by the user input script to make a sketch of the flow domain for sphere-heat-transfer.
-- PJ, 2016-08-31
s = Sketch:new{renderer="svg", projection="xyortho"}
s:set{canvas={0.0,0.0,120.0,120.0}, viewport={Fx-R,Fy-R,Dx+R,Dy+R}}

s:start{file_name="alkandry2014-grid.svg"}
-- s:text{point=Vector3:new{x=-0.007,y=0.022}, text="flow domain for sphere", font_size=16}

s:set{line_width=0.1, fill_colour="green"} -- for flow domain
s:render{surf=psrf0}
s:render{surf=psrf1}

s:set{line_width=0.5}; -- thick line for sphere surface
s:dotlabel{point=p, label="p"}
s:dotlabel{point=q, label="q"}
s:dotlabel{point=a, label="a"}
s:dotlabel{point=b, label="b"}
s:dotlabel{point=c, label="c"}


s:render{path=mid}
s:render{path=ga}
s:render{path=abc0}
s:render{path=gh0}
s:render{path=hc}
s:render{path=abc1}
s:render{path=gh1}

s:finish{}
