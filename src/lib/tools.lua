-- A small set of mathematical tools.
-- Jens Christoph Kunze <j.kunze@uq.edu.au>
--

local tools = {}

function mean(table)
    -- Return arithmetic mean of a table
    table_sum = sum(table)
    table_mean = table_sum / #table

    return table_mean
end

function sum(table)
    -- Return sum of the elements of an array
    local s = 0
    for i,v in ipairs(table) do
        s = s + v
    end
    return s
end

local function round(num)
    -- Return rounded value to the nearest integer
    -- of a given number and delete trailing zeros
    num = math.floor(num + 0.5)
    return num
end

local function round2(num,prec)
    -- Return rounded value to given precision
    -- of a given number and delete trailing zeros
    num = math.floor(num / (10^prec) + 0.5)
    return num * 10^prec
end

local function max(table)
    -- Return max value in a table
    max = table[1]
    max_index = 1
    for i, num in ipairs(table) do
        if table[i] > max then
            max = table[i]
            max_index = i
        end
    end
    return max, max_index
end

local function extreme(table)
    -- Return max absolute value in a table
    extreme = math.abs(table[1])
    extreme_index = 1
    for i, num in ipairs(table) do
        if math.abs(table[i]) > math.abs(extreme) then
            extreme = table[i]
            extreme_index = i
        end
    end
    return extreme, extreme_index
end

local function min(table)
    -- Return min value in a table
    min = table[1]
    min_index = 1
    for i, num in ipairs(table) do
        if table[i] < min then
            min = table[i]
            min_index = i
        end
    end
    return min, min_index
end

local function lin_interp(x_0,x_1,y_0,y_1,x)
    -- Linear interpolation
    local y = y_0 + (y_1 - y_0) / (x_1 - x_0) * (x - x_0)

    return y
end

local function mov_avg(tab,width)
    -- Moving average of a table
    local width_m = math.floor(width / 2)

    -- Copy data for values at either end
    tab_avg = {}
    for i,v in ipairs(tab) do
        tab_avg[i] = v
    end

    if  width_m == width / 2 then
        width_p = width_m - 1
    else
        width_p = width_m
    end

    for i=1+width_m, #tab-width_p do
        local mov_sum = 0
        for j=i-width_m, i+width_p do
            mov_sum = mov_sum + tab[j]
        end
        tab_avg[i] = mov_sum / width
    end

    return tab_avg
end

local function vec_prod(vec1,vec2)
    -- Vector product of two given vectors
    local n = {}
    n[1] = vec1[2] * vec2[3] - vec1[3] * vec2[2]
    n[2] = vec1[3] * vec2[1] - vec1[1] * vec2[3]
    n[3] = vec1[1] * vec2[2] - vec1[2] * vec2[1]

    return n
end

local function factorial(number)
    -- factorial of a given number
    local factorial = 0
    while number > 1 do
        factorial = factorial * number
        number = number - 1
    end
    return factorial
end

local function self_sum(number)
    -- sum equivalent of a factorial
    -- n + n-1 + n-2 + ... + 3 + 2 + 1
    self_sum = math.ceil(number/2) * number + math.fmod(number+1,2) * number/2
    return self_sum
end

local function sleep(seconds)
    -- takes number of seconds as an input and sleeps for the given time
    local sleep_time = os.clock() + seconds
    repeat until os.clock() >= sleep_time
    return 0
end

local function scal_prod(vec1,vec2)
    -- calculates scalar product or dot product of two vectors
    local sp = vec1[1] * vec2[1] + vec1[2] * vec2[2] + vec1[3] * vec2[3]
    return sp
end

local function vec_ang(vec1,vec2)
    -- calculates angle in between two vectors
    local sp = scal_prod(vec1,vec2)
    local l1 = scal_prod(vec1,vec1)
    local l2 = scal_prod(vec2,vec2)
    local ang = math.acos(sp/math.sqrt(l1*l2))
    return ang
end


local function bilin_interp(x, y, val, ip)
    -- bilinear interpolation
    -- x1 < x2, y1 < y2
    -- val1 @ x1y1, val2 @ x2y1, val3 @ x1y2, val4 @ x2y2
    -- ip := interpolation point (x,y)
    local mx1 = (x[2] - ip[1]) / (x[2] - x[1])
    local mx2 = (ip[1] - x[1]) / (x[2] - x[1])
    local my1 = (y[2] - ip[2]) / (y[2] - y[1])
    local my2 = (ip[2] - y[1]) / (y[2] - y[1])
    local res = my1 * (mx1 * val[1] + mx2 * val[2]) + my2 * (mx1 * val[3] + mx2 * val[4])
    return res
end

local function multilin_interp(x, y, val, ip)
    -- multilin interpolation for arbitrary points
    -- weighted average of the given number of points
    -- x= (x1, .. ,xn)
    -- y= (y1, .. ,yn)
    -- val= (val1, .. ,valn)
    -- d1-dn := distance from interpolation point to points 1-n
    local d = {}
    local dt = 0
    for i=1, #x do
        d[#d+1] = math.sqrt((y[i] - ip[2])^2 + (x[i] - ip[1])^2)^(-1)
        dt = dt + d[#d]
    end
    local valp = 0
    for i=1, #x do
        valp = valp + val[i] * d[i] / dt
    end
    return valp
end

local function three_point_plane(p1,p2,p3)
    -- create plane through three points
    -- returns normal vector of plane
    vec1 = {p1[1]-p2[1],p1[2]-p2[2],p1[3]-p2[3]}
    vec2 = {p3[1]-p2[1],p3[2]-p2[2],p3[3]-p2[3]}
    local n = vec_prod(vec1,vec2)
    return n
end

local function rms(table)
    -- root mean square of the elements of an array
    -- square each value first
    local sq = {}
    for i=1, #table do
        sq[i] = table[i]^2
    end

    -- calculate rms
    local res = math.sqrt(sum(sq) / #sq)
    return res
end

local function sort_by_x(t1,t2)
    -- Two tables are sorted according to their x-value
    -- when used with table.sort
    -- Syntax: table.sort(table, sort_by_x)
    local a = t1.x
    local b = t2.x
    return a < b
end

return {mean = mean,
    sum = sum,
    round = round,
    round2 = round2,
    max = max,
    extreme = extreme,
    min = min,
    lin_interp = lin_interp,
    mov_avg = mov_avg,
    vec_prod = vec_prod,
    factorial = factorial,
    self_sum = self_sum,
    sleep = sleep,
    scal_prod = scal_prod,
    vec_ang = vec_ang,
    copy_table = copy_table,
    bilin_interp = bilin_interp,
    multilin_interp = multilin_interp,
    three_point_plane = three_point_plane,
    rms = rms,
    sort_by_x = sort_by_x
    }
