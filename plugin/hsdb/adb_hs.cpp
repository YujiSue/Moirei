#include "../../moirei.h"
///////////////////////////////////////////// slib namespace ////////////////////////////////////////////
using namespace slib;
using namespace slib::sutil;
using namespace slib::sbio;
using namespace slib::sbio::sutil;
/////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////// Data containes ////////////////////////////////////////////
sindex chrindex, refindex;
Array<AnnotInfo> contigs;
Array<GeneInfo> genes, genes2;
Map<String, GeneInfo*> geneIndex;
Map<String, int> geneIndex2;
Array<TranscriptInfo> transcripts;
Map<String, int> trsIndex;
Array<VariantInfo> mutations, variants;
Map<String, int> mutIndex, varIndex;
Array<AnnotInfo> features;
Array<AnnotInfo> proteins;
Map<String, int> protIndex;
Array<MotifInfo> motifs;
Map<int, Array<Pair<String, String>>> crossRef[3];
/////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////// Map to classify ///////////////////////////////////////////
// Transcript type identifier
sindex hs_transcript_type = {
    SI_("mRNA", (int)TRANSCRIPT_TYPE::M_RNA),
    SI_("tRNA", (int)TRANSCRIPT_TYPE::T_RNA),
    SI_("rRNA", (int)TRANSCRIPT_TYPE::R_RNA),
    SI_("ncRNA", (int)TRANSCRIPT_TYPE::NC_RNA),
    SI_("lnc_RNA", (int)TRANSCRIPT_TYPE::LNC_RNA),
    SI_("scRNA", (int)TRANSCRIPT_TYPE::SC_RNA),
    SI_("snRNA", (int)TRANSCRIPT_TYPE::SN_RNA),
    SI_("snoRNA", (int)TRANSCRIPT_TYPE::SNO_RNA),
    SI_("miRNA", (int)TRANSCRIPT_TYPE::MI_RNA),
    SI_("antisense_RNA", (int)TRANSCRIPT_TYPE::AS_RNA),
    SI_("Y_RNA", (int)TRANSCRIPT_TYPE::Y_RNA),
    SI_("RNase_P_RNA", (int)TRANSCRIPT_TYPE::RNASE_P),
    SI_("RNase_MRP_RNA", (int)TRANSCRIPT_TYPE::RNASE_MRP),
    SI_("telomerase_RNA", (int)TRANSCRIPT_TYPE::TEROMERASE_RNA),
    SI_("pseudogenic_transcript", (int)TRANSCRIPT_TYPE::NC_RNA),
    SI_("circular_ncRNA", (int)TRANSCRIPT_TYPE::CIRC_NC_RNA),
    SI_("vault_RNA", (int)TRANSCRIPT_TYPE::VT_RNA),
    SI_("vaultRNA_primary_transcript", (int)TRANSCRIPT_TYPE::VT_RNA)
};
// Structure type identifier
sindex hs_struct_type = {
    SI_("CDS", sbio::CDS),
    SI_("exon", sbio::EXON),
    SI_("miRNA", sbio::PROCESSING)
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////// Func. to show progress ///////////////////////////////////////
inline void _showProgress(GffFile* gff, bool* run) {
    while (!gff->eof() && (*run)) {
        SWrite(slib::DEL * 4, sstr::lfill(S((int)(100 * gff->offset() / gff->size())), ' ', 3), "%");
        sleep(1000);
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/*****************************************************/
/*  Unique function to process each species dataset  */
/*****************************************************/


/*
SDictionary hs_gff = {
    D_("RefSeq", sobj({
        
        D_("biological_region", sannot::FEATURE),
        D_("region", sannot::FEATURE),
        D_("sequence_feature", sannot::FEATURE),
        D_("sequence_alteration_artifact", sannot::FEATURE),
        D_("D_loop", sannot::FEATURE),
        D_("centromere", sannot::FEATURE)
    })),



    D_("Gnomon", sobj({
        
    

        D_("V_gene_segment", sannot::FEATURE)
    })),

    

    D_("Curated Genomic", sobj({

        D_("C_gene_segment", sannot::FEATURE),
        D_("D_gene_segment", sannot::FEATURE),
        D_("J_gene_segment", sannot::FEATURE),
        D_("V_gene_segment", sannot::FEATURE),
        D_("recombination_feature", sannot::FEATURE),
        D_("sequence_feature", sannot::FEATURE),
        D_("enhancer", sannot::FEATURE)
    })),

    D_("RefSeqFE", sobj({
        D_("C_gene_segment", sannot::FEATURE),
        D_("D_gene_segment", sannot::FEATURE),
        D_("J_gene_segment", sannot::FEATURE),
        D_("V_gene_segment", sannot::FEATURE),
        D_("recombination_feature", sannot::FEATURE),
        D_("sequence_feature", sannot::FEATURE),
        D_("enhancer", sannot::FEATURE)
    }))
};
*/
inline void makeRefRenameMap(const SeqList& reference) {
    sfori(reference) {
        refindex[reference[i].attribute["genbank"]] = i;
    }
}
inline void makeGeneIndex(const char* path) {
    SFile f(path);
    SString txt;
    f >> txt;
    auto lines = txt.trim().splitLine();
    lines.removeAt(0);
    genes.resize(lines.size());
    sfor2(lines, genes) {
        auto vals = $_1.split(TAB);
        $_2.geneid = vals[0];
        geneIndex[$_2.geneid] = it.second.ptr();
        $_2.name = vals[1];
        geneIndex[$_2.name] = it.second.ptr();
        $_2.description = vals[2];
        if (vals[3] == "Approved") $_2.type = (int)GENE_TYPE::MISC_GENE;
        else {
            $_2.type = (int)GENE_TYPE::UNAVAILABLE;
            continue;
        }
        if (4 < vals.size() && vals[4].size()) $_2.synonym.append(vals[4].split(","));
        if (5 < vals.size() && vals[5].size()) $_2.synonym.append(vals[5].split(","));
        if (6 < vals.size() && vals[6].size()) $_2.synonym.append(vals[6].split(","));
        if (7 < vals.size() && vals[7].size()) crossRef[0][it.second - genes.begin() + 1].add(Pair<String, String>("genbank", vals[7]));
        if (8 < vals.size() && vals[8].size()) {
            geneIndex[vals[8]] = it.second.ptr();
            crossRef[0][it.second - genes.begin() + 1].add(Pair<String, String>("ncbi", vals[8]));
        }
        if (9 < vals.size() && vals[9].size()) crossRef[0][it.second - genes.begin() + 1].add(Pair<String, String>("ensembl", vals[9]));
    }
}


inline void _setGffData(AnnotInfo* info, GffData* data) {
    info->idx = refindex[data->seqid];
    info->begin = data->begin;
    info->end = data->end;
    info->dir = data->dir;
}
inline void _addGeneData(GffData* data) {
    GeneInfo* info = nullptr;
    auto& symbol = data->attribute["Name"];
    if (geneIndex.hasKey(symbol)) info = geneIndex[symbol]; 
    else {
        geneIndex2.set(symbol, genes2.size());
        genes2.add();
        info = &genes2[-1];
    }
    _setGffData(info, data);
    if (data->type == "pseudogene") info->type = (int)GENE_TYPE::PSEUDO_GENE;
    else if (data->attribute.hasKey("gene_biotype")) {
        auto& type = data->attribute["gene_biotype"];
        if (type.beginWith("protein")) info->type = (int)GENE_TYPE::PROTEIN_CODING;
        else if (type.beginWith("pseudo")) info->type = (int)GENE_TYPE::PSEUDO_GENE;
        else info->type = (int)GENE_TYPE::NON_CODING;
    }
}
inline void _addTranscriptData(GffData* data) {
    if (data->type == "miRNA") return;
    transcripts.add();
    auto& tinfo = transcripts[-1];
    _setGffData(&tinfo, data);
    auto& symbol = data->attribute["gene"];
    if (geneIndex.hasKey(symbol)) {
        tinfo.gene = geneIndex[symbol];
        // Set index as gene record index
        tinfo.idx = tinfo.gene - genes.data() + 1;
    }
    else if (geneIndex2.hasKey(symbol)) {
        tinfo.gene = nullptr;
        // Set index as gene record index
        tinfo.idx = -1;
    }
    else {
        SPrint(symbol, " was not listed.");
        tinfo.gene = nullptr;
        tinfo.idx = -1;
    }
    if (hs_transcript_type.hasKey(data->type))
        tinfo.type = hs_transcript_type[data->type];
    else if (data->type == "primary_transcript") {
        if (data->attribute.hasKey("miRBase"))  tinfo.type = (int)TRANSCRIPT_TYPE::MI_RNA;
        else tinfo.type = (int)TRANSCRIPT_TYPE::NC_RNA;
    }
    else if (data->attribute["pseudo"] == "true")
        tinfo.type = (int)TRANSCRIPT_TYPE::PSEUDO_GENE_TRANSCRIPT;
    else {
        if (data->attribute.hasKey("gbkey") && data->attribute["gbkey"] == "misc_RNA") {}
        //else SPrint(symbol, " was unknown type : ", data->type);
        tinfo.type = (int)TRANSCRIPT_TYPE::MISC_RNA;
    }
    if (data->attribute.hasKey("Name")) 
        tinfo.name = data->attribute["Name"];
    else if (data->attribute.hasKey("ID")) 
        tinfo.name = data->attribute["ID"].replace("rna-", "");
    else {
        //SPrint("No name transcript.", NL, data->attribute);
    }
    trsIndex[tinfo.name] = transcripts.size() - 1;
}
inline void _addStructureData(GffData* data) {
    AnnotInfo info;
    _setGffData(&info, data);
    info.type = hs_struct_type[data->type];
    auto tid = data->attribute["transcript_id"];
    if (tid == "") tid = data->attribute["Parent"].replace("rna-", "");
    if (trsIndex.hasKey(tid)) {
        info.idx = trsIndex[tid] + 1;
        transcripts[info.idx - 1].structures.add(info);
    }
    else SPrint(tid, " was not indexed transcript.");
}

inline void _addFeatureData(GffData* data) {
    AnnotInfo info;
    _setGffData(&info, data);
    if (data->type == "enhancer") info.type = sbio::ENHANCER;
    else if (data->type == "silencer") info.type = sbio::SILENCER;
    else if (data->type == "promoter") info.type = sbio::PROMOTER;
    else if (data->type == "insulator") info.type = sbio::INSULATOR;
    else if (data->type.beginWith("transcription")) info.type = sbio::REGULATOR;
    else if (data->type == "TATA_box") info.type = sbio::TATA_BOX;
    else if (data->type == "protein_binding_site") info.type = sbio::BINDING_SITE;

    else if (data->type.match("repeat")) info.type = sbio::REPEAT_SITE;    

    else if (data->type == "match") {
        SPrint(data->attribute);
    }

    else if (data->type == "biological_region") info.type = sbio::MISC_FEATURE;
    else SPrint(data->type);
    info.name = data->attribute["standard_name"];
    features.add(info);
}


// Load GTF and interpret
void loadGFF(const char *path) {
    SPrint("Loading gff3 data.");
    SWrite("Progress :     ");
    bool run = true;
    GffData* data;
    GffFile gff;
    gff.open(path);
    std::thread th(_showProgress, &gff, &run);
    //
    while ((data = gff.next())) {
        if (!data) break;
        if (data->type.endWith("gene")) _addGeneData(data);
        else if (data->type.endWith("transcript") || data->type.endWith("RNA")) _addTranscriptData(data);
        else if (hs_struct_type.hasKey(data->type)) _addStructureData(data);
        else _addFeatureData(data);
    }
    run = false;
    th.join();
    SPrint("");
}

inline void toVarAnnotInfo(VariantInfo &vinfo, Variant* var, SDictionary& attribute) {
    vinfo.idx = var->pos[0].idx;
    vinfo.begin = var->pos[0].begin;
    vinfo.end = var->pos[0].end;
    vinfo.name = var->varid;
    vinfo.varid = var->varid;
    sfor(var->attribute) {




    }
}
// Load ClinVar data
void loadClinVar(const char* path) {
    sbio::VcfFile vcf(path);
    Variant* var;
    while (var = vcf.next()) {
        variants.add();
        toVarAnnotInfo(variants[-1], var, vcf.header);
    }
}

// Load motifs exported by InterProScan
inline void addMotifs(const char* path) {
    SFile f(path);
    String ln;
    f.readLine(ln);
    while (f) {
        f.readLine(ln);
        if (ln.empty()) continue;
        auto dat = ln.split(TAB);
        motifs.add();
        auto& motif = motifs[-1];
        motif.name = dat[5];
        motif.type = 0;
        if (protIndex.hasKey(dat[0])) motif.idx = protIndex.hasKey(dat[0]);
        else {
            motif.idx = -1;
            SPrint(dat[0], " was not listed protein.");
        }
        motif.begin = dat[6].intValue();
        motif.end = dat[7].intValue();
        motif.motid = dat[4];
        motif.program = dat[3];
    }
}
/*****************************************************/

///////////////////////////////////////////// Main functtion ////////////////////////////////////////////
extern "C" {
    splugin makeDB( SDataBase& db, 
                    const SDictionary& pref, 
                    const SeqList& reference ) {
        try {
            // Make reference index
            refindex = reference.nameIndex();
            // Get additional datasource
            auto attributes = pref["attribute"].parse(";", "=");

            /**************************************************************/
            /*  Unique codes to process each species-specific data set    */
            /**************************************************************/

            // Make a map to convert : Accession IDs => Common symbols
            makeRefRenameMap(reference); 
            

            // Make gene list based on the HGNC IDs list.
            makeGeneIndex(attributes["genelist"]);
            
            // Load GFF data and interpret
            loadGFF(pref["gff"]);

            /*************************************************************/

            ////////// Append records to SQLite DB //////////
            // Make contig table
            {
                SPrint("Write out contig data.");
                auto ctgtbl = db["contig"];
                ctgtbl.prepare().insert();
                sfor(contigs) { ctgtbl.addRecord({ $_.name, $_.idx, $_.begin, $_.end }); }
                ctgtbl.complete();
            }
            // Make gene and synonym table
            {
                SPrint("Write out gene data.");
                auto genetbl = db["gene"];
                genetbl.prepare().insert();
                sfor(genes) {
                    genetbl.addRecord({ snull, $_.name, $_.type, $_.idx, $_.begin, $_.end, $_.dir,
                        $_.geneid, $_.description, sjson::toString($_.attribute) });
                }
                genetbl.complete();
                auto synontbl = db["synonym"];
                synontbl.prepare().insert();
                sfor(genes) {
                    sforeach(syn, $_.synonym) { synontbl.addRecord({ $INDEX(genes) + 1, syn }); }
                }
                synontbl.complete();
            }
            // Make transcript table
            {
                SPrint("Write out transcript data.");
                auto trstbl = db["transcript"];
                trstbl.prepare().insert();
                sfor(transcripts) {
                    trstbl.addRecord({ snull, $_.name, $_.type, $_.idx, $_.begin, $_.end });
                }
                trstbl.complete();
            }
            // Make structure table
            {
                SPrint("Write out structure data.");
                auto strcttbl = db["structure"];
                strcttbl.prepare().insert();
                sfor(transcripts) {
                    sforeach(strct, $_.structures) {
                        strcttbl.addRecord({ strct.type, strct.idx, strct.begin, strct.end });
                    }
                }
                strcttbl.complete();
            }
            // Make mutation table
            {
                SPrint("Write out mutation data.");
                auto muttbl = db["mutation"];
                muttbl.prepare().insert();
                sfor(mutations) {
                    muttbl.addRecord({ snull, $_.name, $_.type, $_.idx, $_.begin, $_.end,
                        $_.varid, sjson::toString($_.attribute) });
                }
                muttbl.complete();
            }
            // Make variant table
            {
                SPrint("Write out variant data.");
                auto vartbl = db["variant"];
                vartbl.prepare().insert();
                sfor(variants) {
                    vartbl.addRecord({ snull, $_.name, $_.type, $_.idx, $_.begin, $_.end,
                        $_.varid, sjson::toString($_.attribute) });
                }
                vartbl.complete();
            }
            // Make feature table
            {
                SPrint("Write out feature data.");
                auto ftrtbl = db["feature"];
                ftrtbl.prepare().insert();
                sfor(features) {
                    ftrtbl.addRecord({ $_.name, $_.type, $_.idx, $_.begin, $_.end, $_.dir });
                }
                ftrtbl.complete();
            }
            // Make protein table
            {
                SPrint("Write out protein data.");
                auto prottbl = db["protein"];
                prottbl.prepare().insert();
                sfor(proteins) {
                    if ($_.idx < 0) continue;
                    $_.end = transcripts[$_.idx - 1].coding().length(true) / 3 - 1;
                    prottbl.addRecord({ snull, $_.name, $_.idx, $_.begin, $_.end });
                }
                prottbl.complete();
            }
            // Make motif table
            {
                SPrint("Write out motif data.");
                if (attributes.hasKey("motif")) addMotifs(attributes["motif"]);
                auto mottbl = db["motif"];
                mottbl.prepare().insert();
                sfor(motifs) {
                    if ($_.idx < 0) continue;
                    mottbl.addRecord({ $_.name, $_.type, $_.idx, $_.begin, $_.end, $_.motid, $_.program });
                }
                mottbl.complete();
            }
            // Make xref table
            {
                SPrint("Write out xref data.");
                auto xreftbl = db["xref"];
                xreftbl.prepare().insert();
                sforin(i, 0, 3) {
                    sfor(crossRef[i]) {
                        sforeach(x, $_.value()) {
                            xreftbl.addRecord({ i + 1, $_.key(), x.first, x.second });
                        }
                    }
                }
                xreftbl.complete();
            }
            /////////////////////////////////////////////////
            return 0;
        }
        catch (Exception ex) {
            ex.print();
            return ex.code;
        }
    }
}
////////////////////////////////////////////////// End //////////////////////////////////////////////////