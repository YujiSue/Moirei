#include "../../moirei.h"

#if defined(WIN_OS) and defined(_WINDLL)
slib::String slib::NL = "\r\n";
#endif

extern "C" {
    splugin getFilterInfo(slib::SDictionary& info) {
        if (!info.hasKey("filter")) info["filter"] = SDictionary();
        info["filter"]["EMS"] = "Frequent substitutions in EMS mutagenesis.";
        info["filter"]["ENU"] = "Frequent substitutions in ENU mutagenesis.";
        return 0;
    }
    splugin customVarFilter(slib::sbio::Variant* var, sobj opts) {
        bool passed = false;
        if (var->flag & slib::sbio::SMALL_VARIANT) {
            if (var->type == SNV) {
                if (opts["mutagen"] == "EMS") {
                    if ((var->attribute["_ref_"] == "G" && var->alt == "A") ||
                        (var->attribute["_ref_"] == "C" && var->alt == "T")) {
                        var->attribute["filter"] = "EMS"; passed = true;
                    }
                }
                if (opts["mutagen"] == "ENU") {
                    if ((var->attribute["_ref_"] == "A" && var->alt == "T") ||
                        (var->attribute["_ref_"] == "A" && var->alt == "G") ||
                        (var->attribute["_ref_"] == "T" && var->alt == "A") ||
                        (var->attribute["_ref_"] == "T" && var->alt == "C") ||
                        (var->attribute["_ref_"] == "G" && var->alt == "A") ||
                        (var->attribute["_ref_"] == "C" && var->alt == "T")) {
                        var->attribute["filter"] = "ENU"; passed = true;
                    }
                }
            }
        }
        if (!passed) var->flag |= NOT_USE_FLAG;
        return 0;
	}
}