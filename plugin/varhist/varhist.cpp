#include "../../src/moirei.h"

#if defined(WIN_OS) and defined(_WINDLL)
slib::String slib::NL = "\r\n";
#endif

extern "C" {
    splugin customVarIO(const char *path, slib::sbio::VarList* vlist, sobj opts) {
        String output = sfs::joinPath(
            sfs::splitPath(path).first,
            sfs::fileName(path, false) + ".hist.json"
        );
        //
        Array<sveci> histdata;
        auto &lgs = vlist->linkageGroups();
        int num = lgs.size(), bin = 1;
        //
        auto bins = opts["bin"].string();
        if (bins[-1] == 'K') {
            bin = bins.clip(0, bins.size() - 1).intValue() * 1000;
        }
        else if (bins[-1] == 'M') {
            bin = bins.clip(0, bins.size() - 1).intValue() * 1000000;
        }
        else bin = bins.intValue();
        //
        histdata.resize(num);
        sfor2(histdata, lgs) {
            $_1.resize(($_2.second - 1) / bin + 1);
        }
        //
        sforeach(var, *vlist) {
            ++histdata[var->pos[0].idx][(var->pos[0].begin-1) / bin];
        }
        //
        sobj histObj = SArray();
        sforin(i, 0, num) {
            sforin(j, 0, histdata[i].size()) {
                histObj.add({
                    D_("block_id", lgs[i].first),
                    D_("start", j * bin + 1),
                    D_("end", (j + 1) * bin),
                    D_("value", histdata[i][j])
                    });
            }
        }
        sjson::save(histObj, output);
        return 0;
	}
}