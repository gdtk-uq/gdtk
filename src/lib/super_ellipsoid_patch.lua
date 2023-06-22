-- ============================================================================
-- Author: Flynn M. Hack
-- Date: 2023-06-22
--
-- The super_ellipsoid_patch is based on the superquadric class of shapes.
-- These shapes are quite important in computer science, and have some really neat 
-- analytical expressions for the surface, moment of inertia and more.
-- As such, these shapes were also used to model nose cones etc in aerospace.
--
-- The patch is formed by scaling points on a cube face so they 
-- lie on the superellipsoid surface. The general surface is defined by
-- the surface function:
--     1 = ((x / ax) ** (2/e2) + (y / ay) ** (2/e2)) ** (e2/e1) + (z / az) ** (2/e1)
--
-- Currently, the following simplified superellipsoid surface function has been implemented:
--     1 = (x) ** (2/e) + (y) ** (2/e) + z ** (2/e)
--
-- As e --> 0, the surface tends to a perfect cube with sharp edges.
--
-- For more info, the "super ellipsoids" wikipedia page is quite good.
-- Otherwise, Faster Calculation of Superquadric shapes by Franklin and Barr 1981   
-- and the textbook Graphics Gems 3 have some nice explanations.
-- 
-- In future, the addition of ax, ay, az can provide some nice scaling along axes.
-- ===========================================================================

super_ellipsoid_patch = {}

function super_ellipsoid_patch:new(o)

    o = o or {}
    setmetatable(o, self)
    self.__index = self
    return o
end

function super_ellipsoid_patch:eval(r,s)

    face = self.face
    a    = self.a
    e    = self.e
    q0   = self.q0
    q1   = self.q1
    q2   = self.q2
    q3   = self.q3

    -- The above makes the assumption of quaternions relative to RHR coord system 
    -- Modelling and Simulation of Aerospace Vehicles (Zipfel 2007).

    if (q0 == nil or q1 == nil or q2 == nil or q3 == nil) then
        q0 = 1 ; q1 = 0.0 ; q2 = 0.0 ; q3 = 0.0 -- (enforce no rotation)
    end 

    -- First, define our unit cube patch of side length 2a.

    if face == "east" then 
        x_c = 1.0; y_c = -1.0 + 2.0*r; z_c = -1.0 + 2.0*s;

    elseif face == "west" then
        x_c = -1.0; y_c = -1.0 + 2.0*r; z_c = -1.0 + 2.0*s;
    
    elseif face == "south" then
        x_c = -1.0 + 2.0*r; y_c = -1.0; z_c = -1.0 + 2.0*s;
    
    elseif face == "north" then
        x_c = -1.0 + 2.0*r; y_c = 1.0; z_c = -1.0 + 2.0*s;
    
    elseif face == "bottom" then
        x_c = -1.0 + 2.0*r; y_c = -1.0 + 2.0*s; z_c = -1.0;
    
    elseif face == "top" then
        x_c = -1.0 + 2.0*r; y_c = -1.0 + 2.0*s; z_c = 1.0;
    
    else
        return nil
    end

    -- Scale our points to ensure they lie on the super-ellipsoid surface.

    f = 1 / ((math.abs(x_c)^(2/e) + math.abs(y_c)^(2/e) + math.abs(z_c)^(2/e))^(e/2))

    x_c = f * x_c
    y_c = f * y_c
    z_c = f * z_c
   
    -- Apply the rotation.
    
    x_cdash = (x_c * (q0*q0 + q1*q1 - q2*q2 - q3*q3) + y_c * (2*(q1*q2 - q0*q3)) + z_c * (2*(q1*q3 + q0*q2)))
    y_cdash = (x_c * (2*(q1*q2 + q0*q3)) + y_c * (q0*q0 - q1*q1 + q2*q2 - q3*q3) + z_c * (2*(q2*q3 - q0*q1)))
    z_cdash = (x_c * (2*(q1*q3 - q0*q2)) + y_c * (2*(q2*q3 + q0*q1)) + z_c * (q0*q0 - q1*q1 - q2*q2 + q3*q3))   

    -- Scale (in future, we could translate each point after scaling).

    x_cdash = x_cdash * a
    y_cdash = y_cdash * a
    z_cdash = z_cdash * a

    return {x = x_cdash, y = y_cdash, z = z_cdash}

end

function super_ellipsoid_patch:create()
    
    _G["mySurf"] = function(r,s)
        return self:eval(r,s)
    end

    patch = LuaFnSurface:new{luaFnName="mySurf"}   
    return patch
end
