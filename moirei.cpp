#include "moirei.h"
using namespace slib;

moir::Param::Param() {}
moir::Param::~Param() {}
void moir::Param::setPref(SDictionary& pref) {

    if (pref["oformat"]) oformat = pref["oformat"];
}

/**
* Moirei implementation
*/
Moirei::Moirei() {}
Moirei::~Moirei() {}
void Moirei::init() {}
Response Moirei::run(const SDictionary& preference) {    
    auto cmd = preference["_cmd_"];
    if (cmd == "convertGenome") return moir::convertGenome(param, preference);
    else if (cmd == "extractGenome") return moir::extractGenome(param, preference);
    else if (cmd == "referenceSummary") return moir::referenceSummary(param, preference);
    else if (cmd == "countGCRatio") return moir::countGCRatio(param, preference);
    else if (cmd == "refSearch") return moir::refSearch(param, preference);
    //else if (cmd == "primerInfo") return moir::primerInfo(preference);
    
    else if (cmd == "geneInfo") return moir::geneInfo(param, preference);
    //else if (cmd == "getTranscript") return moir::getTranscript(preference);
    //else if (cmd == "getProtein") return moir::getProtein(preference);
    else if (cmd == "makeOrthoDB") return moir::makeOrthoDB(param, preference);

    else if (cmd == "gffSummary") return moir::gffSummary(param, preference);
    else if (cmd == "makeAnnotDB") return moir::makeAnnotDB(param, preference);
    //else if (cmd == "bioAnnotation") return moir::bioAnnot(preference);
    //
    else if (cmd == "makeDiseaseDB") return moir::makeDiseaseDB(param, preference);


    //else if (cmd == "getMotif") return moir::getMotif(preference);

    
    else if (cmd == "variantSearch") return moir::varSearch(preference);
    else if (cmd == "variantFilter") return moir::varFilter(param, preference);

    else if (cmd == "templateParam") {
        sfor(preference["_args_"]) {
            if ($_ == "seq") {
                SeqSearchParam ssp;
                sjson::save(ssp.toObj(), sfs::joinPath(preference["outdir"], "seq.param.json"));
            }
            else if ($_ == "align") {
                AlignmentParam ap;
                sjson::save(ap.toObj(), sfs::joinPath(preference["outdir"], "align.param.json"));
            }
            else if ($_ == "var") {
                VarParam vp;
                sjson::save(vp.toObj(), sfs::joinPath(preference["outdir"], "variant.param.json"));
            }
        }
        return Response();
    }

    //else if (cmd == "geneMap") return moir::geneMap(preference);
    //else if (cmd == "transcriptMap") return moir::transcriptMap(preference);
    //else if (cmd == "motiFMap") return moir::motifMap(preference);
    return Response(-1, "Command is not defined.");
}

