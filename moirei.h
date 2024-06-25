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
        String ln; // Line read buffer
        SFile ofile; // Output file object
        IOStream ostream; // Output stream
        String oformat; // Output format
    public:
        Param();
        ~Param();
        void setPref(SDictionary& pref);
    };
    
    //extern Response downloadRef(moir::Param &param, const SDictionary& pref);
    extern Response convertGenome(moir::Param& param, const SDictionary& pref);
    extern Response extractGenome(moir::Param& param, const SDictionary& pref);
    extern Response referenceSummary(moir::Param& param, const SDictionary& pref);
    extern Response countGCRatio(moir::Param& param, const SDictionary& pref);
    //
    extern Response refSearch(moir::Param& param, const SDictionary& pref);
    extern Response primerInfo(const SDictionary& pref);


    extern Response getTranscript(moir::Param& param, const SDictionary& pref);
    
    
    
    //
    extern Response gffSummary(moir::Param& param, const SDictionary& pref);
    extern Response makeAnnotDB(moir::Param& param, const SDictionary& pref);
    extern Response bioAnnot(moir::Param& param, const SDictionary& pref);


    extern Response geneInfo(moir::Param& param, const SDictionary& pref);

    extern Response makeMotifList(const SDictionary& pref);
    
    extern Response getProtein(moir::Param& param, const SDictionary& pref);
    extern Response getMotif(moir::Param& param, const SDictionary& pref);

    extern Response makeOrthoDB(moir::Param& param, const SDictionary& pref);
    extern Response getOrthoGenes(moir::Param& param, const SDictionary& pref);

    extern Response makeDiseaseDB(moir::Param& param, const SDictionary& pref);
    extern Response diseaseRelatedGene(moir::Param& param, const SDictionary& pref);


    //extern Response digestionSite(const SDictionary& pref);

    extern Response varSearch(const SDictionary& pref);
    extern Response varFilter(moir::Param& param, const SDictionary& pref);

    extern Response geneMap(const SDictionary& pref);
    extern Response transcriptMap(const SDictionary& pref);
    extern Response motifMap(const SDictionary& pref);

}
class Moirei {
public:
    moir::Param param;
public:
    Moirei();
    ~Moirei();
    void init();
    Response run(const slib::SDictionary& preference);
};

#endif