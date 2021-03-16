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
    y = y_0 + (y_1 - y_0) / (x_1 - x_0) * (x - x_0)

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
        width_p = width_m + 1
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
    while number > 0 do
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
    self_sum = self_sum}
