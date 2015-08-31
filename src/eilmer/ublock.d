// ublock.d
// Class for unstructured blocks of cells, for use within Eilmer4.
// Peter J. 2014-07-20 first cut.

module ublock;

import std.conv;
import block;

class UBlock: Block {
public:
    this(int id)
    {
	super(id, "dummy UBlock");
    }


    override string toString() const
    {
	char[] repr;
	repr ~= "UBlock(";
	repr ~= "id=" ~ to!string(id);
	repr ~= ")";
	return to!string(repr);
    }

} // end class UBlock
