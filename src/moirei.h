#ifndef MOIREI_H
#define MOIREI_H
#include "sobj.h"
#include "sbioinfo.h"
#include "sapp/scuiapp.h"
using namespace slib;
using namespace slib::smath;
using namespace slib::sio;
using namespace slib::sutil;
using namespace slib::sbio;
using namespace slib::sbio::sio;
using namespace slib::sbio::sutil;
namespace moir {
    class Param {
    public:
        String ln; // Input read buffer
        String outdir, output, oformat;
        stringarray opath;
        SFile ofile; // Output file object
        IOStream ostream; // Output stream
        //
        Response response; // Result container
    public:
        Param();
        ~Param();
        void setPref(const SDictionary& pref);
    };
    extern void exportSeqs(Param& param, SeqList& seqs);
    extern void rearrange(SFigure& fig, float margin);
}
class Moirei {
public:
    moir::Param param;
public:
    Moirei();
    ~Moirei();
    void init();
    void setParam(SDictionary& preference);
    int run(slib::SDictionary& preference);

    Response& convertGenome(const SDictionary& pref);
    Response& makeAnnotDB(const SDictionary& pref);
    Response& makeOrthoDB(const SDictionary& pref);
    Response& makeDiseaseDB(const SDictionary& pref);

    Response& refInfo(const SDictionary& pref);
    Response& gffSummary(const SDictionary& pref);
    Response& geneInfo(const SDictionary& pref);
    Response& orthoGene(const SDictionary& pref);
    Response& diseaseRelated(const SDictionary& pref);
    Response& diseaseInfo(const SDictionary& pref);
    Response& relatedDisease(const SDictionary& pref);
    Response& variantInfo(const SDictionary& pref);

    Response& primerInfo(const SDictionary& pref);
    Response& motifInfo(const SDictionary& pref);
    
    Response& genomeSeq(const SDictionary& pref);
    Response& geneSeq(const SDictionary& pref);
    Response& transcriptSeq(const SDictionary& pref);
    Response& getComplement(const SDictionary& pref);
    Response& getTranslated(const SDictionary& pref);

    Response& refSearch(const SDictionary& pref);
    Response& varSearch(const SDictionary& pref);
    Response& bioAnnot(const SDictionary& pref);

    Response& countGCRatio(const SDictionary& pref);
    Response& varFilter(const SDictionary& pref);

    Response& geneMap(const SDictionary& pref);
    Response& motifMap(const SDictionary& pref);
};

#endif