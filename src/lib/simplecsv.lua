--[[
https://nocurve.com/2014/03/05/simple-csv-read-and-write-using-lua/
amnon.david@gmail.com

-------------------------------
Save this file as simplecsv.lua
-------------------------------

Example use: Suppose file csv1.txt is:

1.23,70,hello
there,9.81,102
x,y,,z
8,1.243,test

Then the following
-------------------------------------
local csvfile = require "simplecsv"
local m = csvfile.read('./csv1.txt') -- read file csv1.txt to matrix m
print(m[2][3])                       -- display element in row 2 column 3 (102)
m[1][3] = 'changed'                  -- change element in row 1 column 3
m[2][3] = 123.45                     -- change element in row 2 column 3
csvfile.write('./csv2.txt', m)       -- write matrix to file csv2.txt
-------------------------------------

will produce file csv2.txt with the contents:

1.23,70,changed
there,9.81,123.45
x,y,,z
8,1.243,test

the read method takes 4 parameters:
path: the path of the CSV file to read - mandatory
sep: the separator character of the fields. Optionsl, defaults to ','
tonum: whether to convert fields to numbers if possible. Optional. Defaults to true
null: what value should null fields get. Optional. defaults to ''
]]

local scsv = {}

---------------------------------------------------------------------
function string:split(sSeparator, nMax, bRegexp)
    if sSeparator == '' then
        sSeparator = ','
    end

    if nMax and nMax < 1 then
        nMax = nil
    end

    local aRecord = {}

    if self:len() > 0 then
        local bPlain = not bRegexp
        nMax = nMax or -1

        local nField, nStart = 1, 1
        local nFirst,nLast = self:find(sSeparator, nStart, bPlain)
        while nFirst and nMax ~= 0 do
            aRecord[nField] = self:sub(nStart, nFirst-1)
            nField = nField+1
            nStart = nLast+1
            nFirst,nLast = self:find(sSeparator, nStart, bPlain)
            nMax = nMax-1
        end
        aRecord[nField] = self:sub(nStart)
    end

    return aRecord
end

---------------------------------------------------------------------
function scsv.read(path, sep, tonum, null)
    tonum = tonum or true
    sep = sep or ','
    null = null or ''
    local csvFile = {}
    local file = assert(io.open(path, "r"))
    for line in file:lines() do
        fields = line:split(sep)
        if tonum then -- convert numeric fields to numbers
            for i=1,#fields do
                local field = fields[i]
                if field == '' then
                    field = null
                end
                fields[i] = tonumber(field) or field
            end
        end
        table.insert(csvFile, fields)
    end
    file:close()
    return csvFile
end

---------------------------------------------------------------------
function scsv.write(path, data, sep)
    sep = sep or ','
    local file = assert(io.open(path, "w"))
    for i=1,#data do
        for j=1,#data[i] do
            if j>1 then file:write(sep) end
            file:write(data[i][j])
        end
        file:write('\n')
    end
    file:close()
end

---------------------------------------------------------------------

return scsv
