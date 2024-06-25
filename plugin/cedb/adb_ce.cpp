#include "../../moirei.h"

using namespace slib;
using namespace slib::sbio;

sindex chrindex;
Array<AnnotInfo> contigs;
Array<GeneInfo> genes;
Map<String, GeneInfo*> geneIndex;
Array<TranscriptInfo> transcripts;
Map<String, int> trsIndex;
Array<VariantInfo> mutations, variants;
Map<String, int> mutIndex, varIndex;
Array<AnnotInfo> features;
Array<AnnotInfo> proteins;
Map<String, int> protIndex;
Array<MotifInfo> motifs;
Map<int, Array<Pair<String, String>>> crossRef[3];

sindex ce_transcript_type = {
    S_i("mRNA", (int)TRANSCRIPT_TYPE::M_RNA),
    S_i("tRNA", (int)TRANSCRIPT_TYPE::T_RNA),
    S_i("rRNA", (int)TRANSCRIPT_TYPE::R_RNA),
    S_i("ncRNA", (int)TRANSCRIPT_TYPE::NC_RNA),
    S_i("circular_ncRNA", (int)TRANSCRIPT_TYPE::CIRC_NC_RNA),
    S_i("nc_primary_transcript", (int)TRANSCRIPT_TYPE::NC_RNA),
    S_i("lincRNA", (int)TRANSCRIPT_TYPE::LINC_RNA),
    S_i("miRNA_primary_transcript", (int)TRANSCRIPT_TYPE::NC_RNA),
    S_i("pre_miRNA", (int)TRANSCRIPT_TYPE::NC_RNA),
    S_i("miRNA", (int)TRANSCRIPT_TYPE::MI_RNA),
    S_i("antisense_RNA", (int)TRANSCRIPT_TYPE::AS_RNA),
    S_i("piRNA", (int)TRANSCRIPT_TYPE::PI_RNA),
    S_i("snRNA", (int)TRANSCRIPT_TYPE::SN_RNA),
    S_i("snoRNA", (int)TRANSCRIPT_TYPE::SNO_RNA),
    S_i("scRNA", (int)TRANSCRIPT_TYPE::SC_RNA),
    S_i("pseudogenic_transcript", (int)TRANSCRIPT_TYPE::NC_RNA),
    S_i("pseudogenic_tRNA", (int)TRANSCRIPT_TYPE::NC_RNA),
    S_i("pseudogenic_rRNA", (int)TRANSCRIPT_TYPE::NC_RNA)
};
sindex ce_struct_type = {
    S_i("CDS", sbio::CDS),
    S_i("five_prime_UTR", sbio::UTR5),
    S_i("three_prime_UTR", sbio::UTR3),
    S_i("exon", sbio::EXON)
};
//
inline void makeGeneIndex(const char* path) {
    SFile f(path);
    SString txt;
    f >> txt;
    auto lines = txt.splitLine();
    genes.resize(lines.size());
    stringarray data;
    sfor2(lines, genes) {
        data = $_1.split(",");
        if (data.size() < 2 || !(data[1].beginWith("WBG"))) continue;
        geneIndex[data[1]] = &($_2);
        $_2.geneid = data[1];
        $_2.name = data[2].empty() ? data[3] : data[2];
        geneIndex[$_2.name] = &($_2);
        geneIndex[data[3]] = &($_2);
        $_2.synonym = data[3].split(",");
        $_2.type = (suint)(data[4] == "Dead" ? GENE_TYPE::UNAVAILABLE : GENE_TYPE::MISC_GENE);
    }
}
inline void addDescription(const char* path) {
    SFile f(path);
    String ln;
    while (f) {
        f.readLine(ln, true);
        if (ln.beginWith("#")) continue;
        auto dat = ln.split(TAB);
        if (dat[0].empty() || !geneIndex.hasKey(dat[0])) continue;
        auto& gene = *geneIndex[dat[0]];
        // 
        while (f) {
            f.readLine(ln, true);
            if (ln.empty()) continue;
            else if (ln[0] == '=') break;
            gene.description << ln;
        }
        gene.description.trim();
    }
}
inline void setXrefGenes(const char* path) {
    SFile f(path);
    String ln;
    f.readLine(ln);
    while (f) {
        f.readLine(ln, true);
        if (ln.beginWith("#")) continue;
        auto dat = ln.split(TAB);
        if (dat.size() < 6 || dat[4] != "live" || !geneIndex.hasKey(dat[5])) continue;
        auto gene = geneIndex[dat[5]];
        crossRef[0][(gene - genes.data() + 1)].add(Pair<String, String>("ncbi", dat[2]));
    }
}
inline void _setGffData(AnnotInfo* info, GffData* data) {
    info->idx = chrindex[data->seqid];
    info->begin = data->begin;
    info->end = data->end;
    info->dir = data->dir;
}
inline void _addCtgData(GffData* data) {
    contigs.add();
    _setGffData(&contigs[-1], data);
    contigs[-1].name = data->attribute["Name"];
}
inline void _addGeneData(GffData* data) {
    auto& gid = data->attribute["ID"];
    gid.clip(gid.find(":") + 1);
    if (!geneIndex.hasKey(gid)) {
        SPrint(gid, " was not listed.");
        return;
    }
    auto info = geneIndex[gid];
    _setGffData(info, data);
    if (data->source.match("transposon")) info->type = (int)GENE_TYPE::TRANSPOSON;
    else {
        auto& type = data->attribute["biotype"];
        if (type.beginWith("protein")) info->type = (int)GENE_TYPE::PROTEIN_CODING;
        else if (type.beginWith("pseudo")) info->type = (int)GENE_TYPE::PSEUDO_GENE;
        else info->type = (int)GENE_TYPE::NON_CODING;
    }
}
inline void _addTranscriptData(GffData* data) {
    transcripts.add();
    auto& tinfo = transcripts[-1];
    _setGffData(&tinfo, data);
    auto& gid = data->attribute["Parent"];
    gid.clip(gid.find(":") + 1);
    if (geneIndex.hasKey(gid)) {
        tinfo.gene = geneIndex[gid];
        // Set index as gene record index
        tinfo.idx = tinfo.gene - genes.data() + 1;
    }
    else {
        SPrint(gid, " was not listed.");
        tinfo.gene = nullptr;
        tinfo.idx = -1;
    }
    if (ce_transcript_type.hasKey(data->type))
        tinfo.type = ce_transcript_type[data->type];
    else tinfo.type = (int)TRANSCRIPT_TYPE::MISC_RNA;
    tinfo.name = data->attribute["Name"];
    trsIndex[tinfo.name] = transcripts.size() - 1;
    //
    if (tinfo.type == (int)TRANSCRIPT_TYPE::M_RNA && data->attribute.hasKey("wormpep")) {
        proteins.add();
        auto& proInfo = proteins[-1];
        // Set index as transcript record index
        proInfo.idx = transcripts.size();
        proInfo.name = data->attribute["wormpep"];
        proInfo.begin = 1;
        proInfo.end = -1;
        if (data->attribute.hasKey("uniprot_id"))
            crossRef[2][proteins.size()].add("uniprot", data->attribute["uniprot_id"]);
        auto tnames = tinfo.name.split(".");
        if (tnames[-1].match(REG("/\\d+/"))) tnames.resize(tnames.size() - 1);
        protIndex[toString(tnames, ".")] = proteins.size();
    }
}
inline void _addStructData(GffData* data) {
    if (!ce_struct_type.hasKey(data->type)) return;
    AnnotInfo info;
    _setGffData(&info, data);
    info.type = ce_struct_type[data->type];
    auto tids = data->attribute["Parent"].split(",");
    sfor(tids) {
        $_.clip($_.find(":") + 1);
        // Set index as transcript record index
        if (trsIndex.hasKey($_)) {
            info.idx = trsIndex[$_] + 1;
            transcripts[info.idx - 1].structures.add(info);
        }
        else SPrint($_, " was not indexed transcript.");
    }
}
inline void _setVarType(VariantInfo* info, GffData* data) {
    if (data->type == "SNP" || data->type == "substitution" || 
        data->type.match("alter") || data->type.beginWith("point")) {
        info->type = SNV;
        if (data->attribute.hasKey("substitution")) {
            auto& sub = data->attribute["substitution"];
            info->attribute["alt"] = sub.substring(sub.find("/") + 1);
            if (1 < info->attribute["alt"].size()) info->type = MNV;
        }
        if (data->attribute.hasKey("aachange")) {
            info->attribute["substitution"] = data->attribute["aachange"];
        }
    }
    else if (data->type.beginWith("del")) info->type = DELETION;
    else if (data->type.beginWith("ins")) {
        info->type = INSERTION;
        if (data->attribute.hasKey("insertion")) {
            info->attribute["alt"] = data->attribute["insertion"];
        }
    }
    else if (data->type.beginWith("complex")) {
        if (data->attribute.hasKey("insertion")) {
            info->attribute["alt"] = data->attribute["insertion"];
            if (info->length(true) == info->attribute["alt"].size()) info->type = MNV;
            else info->type = DELETION | INSERTION;
        }
        else info->type = DELETION | INSERTION;
    }
    else if (data->type.match("dup")) info->type = DUPLICATION;
}
inline void _setVarConseq(VariantInfo* info, sattribute& attr) {
    if (!attr.hasKey("consequence")) return;
    auto conseq = sstr::toLower(attr["consequence"]);
    auto type = info->type & 0xFF;
    if (conseq.match("missense")) info->type |= ((MISSENSE << 16) | (CDS << 8));
    else if (conseq.match("start")) {
        if (type == DELETION) info->type |= ((HEAD_LESION << 16) | ((CDS | UTR5) << 8));
        else info->type |= ((MISSENSE << 16) | (CDS << 8));
    }
    else if (conseq.match("stop")) {
        if (type == DELETION) info->type |= ((TAIL_LESION << 16) | ((CDS | UTR3) << 8));
        else info->type |= ((MISSENSE << 16) | (CDS << 8));
    }
    else if (conseq.match("synonymous")) info->type |= ((NONSENSE << 16) | (CDS << 8));
    else if (conseq.match("exon_variant")) {
        if (type == SNV || type == MNV) info->type |= (SUBSTITUTION << 16); else info->type |= (INDEL << 16);
        info->type |= (EXON << 8);
    }
    else if (conseq.match("5_prime")) {
        if (type == SNV || type == MNV) info->type |= (SUBSTITUTION << 16); else info->type |= (INDEL << 16);
        info->type |= (UTR5 << 8);
    }
    else if (conseq.match("3_prime")) {
        if (type == SNV || type == MNV) info->type |= (SUBSTITUTION << 16); else info->type |= (INDEL << 16);
        info->type |= (UTR3 << 8);
    }
    else if (conseq.match("splice")) {
        if (type == SNV || type == MNV) info->type |= (SUBSTITUTION << 16); else info->type |= (INDEL << 16);
        info->type |= (SPLICE_SITE << 8);
    }
    else if (conseq.match("intron")) {
        if (type == SNV || type == MNV) info->type |= (SUBSTITUTION << 16); else info->type |= (INDEL << 16);
        info->type |= (INTRON << 8);
    }
    else if (conseq.match("coding_sequence_variant") || conseq.match("protein_altering")) {
        if (type == SNV || type == MNV) info->type |= (SUBSTITUTION << 16); else info->type |= (INDEL << 16);
        info->type |= (CDS << 8);
    }
    else if (conseq.match("ablation")) info->type |= ((NULL_MUT << 16) | (0xFF << 8));
    else if (conseq.match("frameshift")) info->type |= ((FRAME_SHIFT << 16) | (CDS << 8));
    else if (conseq.match("inframe")) info->type |= ((IN_FRAME << 16) | (CDS << 8));
}
inline void _setVarAttr(VariantInfo* info, sattribute& attr) {
    if (attr.hasKey("strain")) info->attribute["strain"] = attr["strain"];

    if (attr.hasKey("vep_impact")) 
        info->attribute["effect"]["evep"] = attr["vep_impact"];
    if (attr.hasKey("polyphen")) 
        info->attribute["effect"]["pph"] = attr["polyphen"];
    if (attr.hasKey("sift")) 
        info->attribute["effect"]["sift"] = attr["sift"];
    if (attr.hasKey("production_method")) 
        info->attribute["method"] = attr["production_method"];
}
inline void _setVarData(VariantInfo *info, GffData* data) {
    _setGffData(info, data);
    info->name = data->attribute["public_name"];
    info->varid = data->attribute["variation"];
    _setVarType(info, data);
    _setVarConseq(info, data->attribute);
    _setVarAttr(info, data->attribute);
}
inline void _addMutData(GffData* data) {
    mutations.add();
    auto& minfo = mutations[-1];
    _setVarData(&minfo, data);
    mutIndex[minfo.name] = mutations.size() - 1;
}
inline void _addVarData(GffData* data) {
    variants.add();
    auto& vinfo = variants[-1];
    _setVarData(&vinfo, data);    
    varIndex[vinfo.varid] = variants.size() - 1;
}
inline void _addFtrData(GffData* data, suint t, const char *s) {
    features.add();
    auto& finfo = features[-1];
    _setGffData(&finfo, data);
    finfo.type = t;
    finfo.name = s;
}

inline void _showProgress(GffFile* gff, bool *run) {
    while (!gff->eof() && (*run)) {
        SWrite(slib::DEL * 4, sstr::lfill(S((int)(100 * gff->offset() / gff->size())), ' ', 3), "%");
        sleep(1000);
    }
}
inline void loadGFF(const char* path) {
    SPrint("Loading gff3 data.");
    SWrite("Progress :     ");
    bool run = true;
    GffData* data;
    GffFile gff;
    gff.open(path);
    std::thread th(_showProgress, &gff, &run);
    //
    while (data = gff.next()) {
        if (!data) break;
        if (data->source == "Genomic_canonical" &&
            data->type == "assembly_component") _addCtgData(data);
        else if (data->source.beginWith("WormBase")) {
            if (data->type == "gene") _addGeneData(data);
            else if (data->type.match("RNA") || data->type.match("transcript")) _addTranscriptData(data);
            else _addStructData(data);
        }
        else if (data->source.endWith("_pmap_position")) {
            auto& note = data->attribute["Note"];
            auto& gid = data->attribute["ID"];
            gid.clip(gid.find(":") + 1);
            if (geneIndex.hasKey(gid) && note.match("cM"))
                geneIndex[gid]->attribute["gmap"] = note.substring(0, note.find("cM")).trim();
        }
        else if (data->source == "landmark") {
            auto& gid = data->attribute["ID"];
            gid.clip(gid.find(":") + 1);
            geneIndex[gid]->attribute["landmark"] = true;
        }
        else if (data->source.match("Polymorphism")) _addVarData(data);
        else if (data->source == "Allele" ||
            data->source == "CGH_allele" ||
            data->source == "KO_consortium" ||
            data->source == "Mos_insertion_allele" ||
            data->source == "NBP_knockout" ||
            data->source == "Variation_project" ||
            data->source == "Million_mutation") _addMutData(data);
        else if (data->source.beginWith("pCoF")) {
            auto& vid = data->attribute["variation"];
            if (mutIndex.hasKey(vid)) mutations[mutIndex[vid]].attribute["phenotype"] = true;
            if (varIndex.hasKey(vid)) variants[varIndex[vid]].attribute["phenotype"] = true;
        }
        else if (data->source.beginWith("Balanced")) {
            auto& name = data->attribute["balancer"];
            name.clip(name.find(":") + 1);
            _addFtrData(data, BALANCED_SITE, name);
        }
        else if (data->source.beginWith("binding_site") ||
            data->source.beginWith("TF_binding")) {
            auto& name = data->attribute.hasKey("tf_name") ? data->attribute["tf_name"] : data->attribute["Name"];
            _addFtrData(data, BINDING_SITE, name);
        }
        else if (data->source == "histone_binding_site_region")
            _addFtrData(data, HISTONE_SITE, data->attribute["Name"]);
        else if (data->source == "miRanda")
            _addFtrData(data, MI_RNA_BINDING, data->attribute["Note"].split(" ")[-1]);
        else if (data->source == "enhancer")
            _addFtrData(data, ENHANCER, data->attribute["Name"]);
        else if (data->source == "operon") 
            _addFtrData(data, OPERON, data->attribute["Name"]);
        else if (data->source == "promoter") 
            _addFtrData(data, PROMOTER, data->attribute["Name"]);
        else if (data->source == "regulatory_region") 
            _addFtrData(data, REGULATOR, data->attribute["Name"]);
        else if (data->source == "TSS_region") _addFtrData(data, TSS_SITE, data->attribute["Name"]);
        else if (data->source.match(REG("/SL[12]/"))) _addFtrData(data, SPLICE_LEADER, data->source);
    }
    run = false;
    th.join();
    SPrint("");
}
inline void addMT(const SeqList& reference) {
    contigs.add();
    auto& info = contigs[-1];
    info.idx = (int)reference.size() - 1;
    info.begin = 1;
    info.end = (int)reference[-1].length();
    info.name = "MTCE";
}
inline void addNBPMutant(const char* path) {
    SFile f(path);
    String ln;
    f.readLine(ln);
    VariantInfo *info;
    while (f) {
        f.readLine(ln);
        if (ln.empty()) continue;
        auto dat = ln.split(",");
        if (mutIndex.hasKey(dat[0])) info = &mutations[mutIndex[dat[0]]];
        else {
            mutations.add();
            info = &mutations[-1];
        }
        info->name = dat[0];
        info->idx = chrindex[dat[1]];
        info->begin = dat[2].intValue();
        info->end = dat[3].intValue();
        info->attribute["alt"] = dat[4];
        //info->type = dat[5].uintValue();
    }
}
inline void addTMBalacer(const char *path) {
    SFile f(path);
    String ln;
    while (f) {
        f.readLine(ln);
        if (ln.empty()) continue;
        auto dat = ln.split(TAB);
        features.add();
        auto& ftr = features[-1];
        ftr.type = sbio::BALANCED_SITE;
        ftr.name = dat[1];
        ftr.idx = dat[2].intValue();
        ftr.begin = dat[3].intValue();
        ftr.end = dat[4].intValue();
    }
}
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
        if (protIndex.hasKey(dat[0])) motif.idx = protIndex[dat[0]];
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
extern "C" {
    splugin makeDB(SDataBase& db, const SDictionary& pref, const SeqList& reference) {
        try {
            chrindex = reference.nameIndex();
            makeGeneIndex(pref["gene-list"]);
            //
            auto suppl = pref.hasKey("supplementary") ? pref["supplementary"].parse(";", "=") : sattribute();
            if (suppl.hasKey("description")) addDescription(suppl["description"]);
            if (suppl.hasKey("xgenes")) setXrefGenes(suppl["xgenes"]);
            //
            loadGFF(pref["gff"]);

            // Make contig table
            {
                SPrint("Write out contig data.");
                addMT(reference);
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
                if (suppl.hasKey("tm")) addNBPMutant(suppl["tm"]);
                auto muttbl = db["mutation"];
                muttbl.prepare().insert();
                sfor(mutations) {
                    muttbl.addRecord({ $_.name, $_.type, $_.idx, $_.begin, $_.end,
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
                    vartbl.addRecord({ $_.name, $_.type, $_.idx, $_.begin, $_.end, 
                        $_.varid, sjson::toString($_.attribute) });
                }
                vartbl.complete();
            }
            // Make feature table
            {
                SPrint("Write out feature data.");
                if (suppl.hasKey("balancer")) addTMBalacer(suppl["balancer"]);
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
                if (suppl.hasKey("motif")) addMotifs(suppl["motif"]);
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
            return 0;
        }
        catch (Exception ex) { 
            ex.print(); 
            return ex.code; 
        }
	}
}