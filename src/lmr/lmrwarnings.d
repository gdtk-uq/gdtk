/**
 * Module to collect warnings in one place.
 *
 * Author: RJG
 * Date: 2024-07-18
 */

module lmr.lmrwarnings;

import std.format;

enum LmrWarning {
    useOfExtremaClipping,
}

string warningMessage(LmrWarning wrn)
{
    final switch (wrn) {
    case LmrWarning.useOfExtremaClipping:
        return
format(`
--- WARNING %02d ---
'extrema_clipping' is set to true.
This is NOT recommended in steady_state mode.
Suggestion: set 'config.extrema_clipping = false' in input script.
------------------
`, LmrWarning.useOfExtremaClipping);
    }
}
