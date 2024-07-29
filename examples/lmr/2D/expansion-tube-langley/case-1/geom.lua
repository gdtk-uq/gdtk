-- geom.lua
--
-- Geometric parameters describing three segments of the main tube.
-- Units are metres.
L_d = 2.44;  R_d = 0.08255 -- driver tube
L_i = 7.49;  R_i = 0.0762  -- intermediate tube
L_a = 14.13; R_a = R_i     -- acceleration tube
L_total = L_d + L_i + L_a
--
-- Detailed geometry as quad patches.
A0 = Vector3:new{x=0.0, y=0.0}
A1 = Vector3:new{x=0.0, y=R_i}
A2 = Vector3:new{x=0.0, y=R_d}
B0 = Vector3:new{x=L_d, y=0.0}
B1 = Vector3:new{x=L_d, y=R_i}
B2 = Vector3:new{x=L_d, y=R_d}
C0 = Vector3:new{x=L_d+L_i, y=0.0}
C1 = Vector3:new{x=L_d+L_i, y=R_i}
D0 = Vector3:new{x=L_d+L_i+L_a, y=0.0}
D1 = Vector3:new{x=L_d+L_i+L_a, y=R_a}
