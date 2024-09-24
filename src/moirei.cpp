#include "moirei.h"
moir::Param::Param() {}
moir::Param::~Param() {}
void moir::Param::setPref(const SDictionary& pref) {
    response.clear();
    if (pref.hasKey("outdir")) outdir = pref["outdir"];
    if (pref.hasKey("oformat")) oformat = pref["oformat"];
    if (pref.hasKey("output")) {
        output = sfs::absolutePath(pref["output"]);
        if (outdir.empty()) outdir = sfs::isDir(output) ? output : sfs::splitPath(output).first;
        if (oformat.empty() || oformat == "auto") oformat = sfs::extension(output);
    }
}
void moir::exportSeqs(moir::Param& param, SeqList& seqs) {
    if (seqs.empty()) return;
    if (param.oformat == "_obj_") {
        sobj seqobj = SArray();
        sfor(seqs) seqobj.add($_.toObj());
        param.response.attribute["result"] = seqobj;
    }
    else if (param.oformat == "fa") {
        if (param.output.empty()) param.ostream.setStrOStream(param.response.output);
        else {
            if (sfs::extension(param.output) != param.oformat) param.output += "." + param.oformat;
            SPrint("Save to '", param.output, "'.");
            param.response.output = param.output;
            param.ofile.open(param.output, MAKE);
            param.ostream.setFileOStream(param.ofile);
        }
        writeFa(param.ostream, seqs);
    }
    else {
        if (seqs.size() == 1) {
            if (param.output.empty()) param.ostream.setStrOStream(param.response.output);
            else {
                if (sfs::extension(param.output) != param.oformat) param.output += "." + param.oformat;
                SPrint("Save to '", param.output, "'.");
                param.response.output = param.output;
                param.ofile.open(param.output, MAKE);
                param.ostream.setFileOStream(param.ofile);
            }
            sbio::sio::writeSeq(param.ostream, seqs[0], param.oformat);
        }
        else {
            if (param.outdir.empty()) param.ostream.setStrOStream(param.response.output);
            sfori(seqs) {
                if (param.outdir.size()) {
                    param.ofile.open(sfs::joinPath(param.outdir, "-" + S(i + 1) + "." + param.oformat), MAKE);
                    SPrint("Save to '", param.ofile.path(), "'.");
                    param.ostream.setFileOStream(param.ofile);
                }
                sbio::sio::writeSeq(param.ostream, seqs[i], param.oformat);
                if (param.ofile.isOpened()) param.ofile.close();
            }
        }
    }
}

void moir::rearrange(SFigure& fig, float margin) {
    if (fig.count() < 2) return;
    while (true) {
        bool overlap = false;
        sforin(it, fig.begin(), fig.end() - 1) {
            auto nxt = it + 1;
            while (nxt < fig.end()) {
                if ($_.overlap(*nxt)) {
                    (*nxt).translate(slib::svec2f(0.f, $_.boundary().height + margin));
                    overlap = true;
                }
                ++nxt;
            }
        }
        if (!overlap) break;
    }
}


/**
* Moirei implementation
*/
Moirei::Moirei() {}
Moirei::~Moirei() {}
void Moirei::init() {}
int Moirei::run(SDictionary& preference) {
    auto cmd = preference["_cmd_"];
    if (cmd == "ConvertGenome") return convertGenome(preference).code;
    else if (cmd == "MakeAnnotDB") return makeAnnotDB(preference).code;
    else if (cmd == "MakeOrthoDB") return makeOrthoDB(preference).code;
    else if (cmd == "MakeDiseaseDB") return makeDiseaseDB(preference).code;

    else if (cmd == "RefInfo") return refInfo(preference).code;
    else if (cmd == "GffSummary") return gffSummary(preference).code;
    else if (cmd == "GeneInfo") return geneInfo(preference).code;
    else if (cmd == "OrthoGene") return orthoGene(preference).code;
    else if (cmd == "DiseaseRelated") return diseaseRelated(preference).code;

    else if (cmd == "DiseaseInfo") return diseaseInfo(preference).code;
    else if (cmd == "RelatedDisease") return relatedDisease(preference).code;


    else if (cmd == "PrimerInfo") return primerInfo(preference).code;

    else if (cmd == "GenomeSeq") genomeSeq(preference).code;
    else if (cmd == "GeneSeq") return geneSeq(preference).code;
    else if (cmd == "TranscriptSeq") return transcriptSeq(preference).code;
    else if (cmd == "GetComplement") return getComplement(preference).code;
    else if (cmd == "GetTranslated") return getTranslated(preference).code;

    else if (cmd == "RefSearch") return refSearch(preference).code;
    else if (cmd == "VariantSearch") return varSearch(preference).code;
    else if (cmd == "BioAnnotation") return bioAnnot(preference).code;

    else if (cmd == "countGCRatio") return countGCRatio(preference).code;
    else if (cmd == "variantFilter") return varFilter(preference).code;
    
    else if (cmd == "GeneMap") return geneMap(preference).code;
    else if (cmd == "MotifMap") return motifMap(preference).code;

    return 0;
}

