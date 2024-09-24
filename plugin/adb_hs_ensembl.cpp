#ifdef HS_E
#include "moirei.h"
using namespace slib;
using namespace slib::sutil;
using namespace slib::sbio;
using namespace slib::sbio::sannot;

String ln;
Array<AnnotInfo> features;
Array<GeneInfo> genes;
Array<TranscriptInfo> transcripts;
Array<VariantInfo> mutations, variants;

Map<String, GeneInfo*> gene_index;
Map<String, sinteger> transcript_index;

SDictionary hs_gff = {
    D_("ensembl", sobj({
        D_("gene", GENE),
        D_("ncRNA_gene", GENE),
        D_("pseudogene", GENE),
        
        D_("mRNA", TRANSCRIPT),
        D_("tRNA", TRANSCRIPT),
        D_("rRNA", TRANSCRIPT),
        D_("ncRNA", TRANSCRIPT),
        D_("lnc_RNA", TRANSCRIPT),
        D_("scRNA", TRANSCRIPT),
        D_("snRNA", "transcript"),
        D_("snoRNA", "transcript"),
        D_("pseudogenic_transcript", "transcript"),
        D_("transcript", "transcript"),
        D_("unconfirmed_transcript", "transcript"),
        
        D_("CDS", "structure"),
        D_("five_prime_UTR", "structure"),
        D_("three_prime_UTR", "structure"),
        D_("exon", "structure")
    })),
    D_("ensembl_havana", sobj({
        D_("gene", "gene"),
        D_("ncRNA_gene", "gene"),
        D_("bidirectional_promoter_lncRNA", "gene"),
        D_("pseudogene", "gene"),
        
        D_("mRNA", TRANSCRIPT),
        D_("lnc_RNA", TRANSCRIPT),
        D_("C_gene_segment", TRANSCRIPT),
        D_("V_gene_segment", "transcript"),
        D_("three_prime_overlapping_ncrna", "transcript"),
        D_("pseudogenic_transcript", "transcript"),
        D_("transcript", "transcript"),
        D_("unconfirmed_transcript", "transcript"),
        
        D_("CDS", "structure"),
        D_("five_prime_UTR", "structure"),
        D_("three_prime_UTR", "structure"),
        D_("exon", "structure")
    })),
    
    D_("havana", sobj({
        D_("gene", "gene"),
        D_("ncRNA_gene", "gene"),
        D_("bidirectional_promoter_lncRNA", "gene"),
        D_("pseudogene", "gene"),
        
        D_("mRNA", "transcript"),
        D_("lnc_RNA", "transcript"),
        D_("scRNA", "transcript"),
        D_("snRNA", "transcript"),
        D_("snoRNA", "transcript"),
        D_("C_gene_segment", "transcript"),
        D_("D_gene_segment", "transcript"),
        D_("J_gene_segment", "transcript"),
        D_("V_gene_segment", "transcript"),
        D_("three_prime_overlapping_ncrna", "transcript"),
        D_("pseudogenic_transcript", "transcript"),
        D_("transcript", "transcript"),
        D_("unconfirmed_transcript", "transcript"),
        D_("vaultRNA_primary_transcript", "transcript"),
        
        D_("CDS", "structure"),
        D_("five_prime_UTR", "structure"),
        D_("three_prime_UTR", "structure"),
        D_("exon", "structure")
    })),
    D_("insdc", sobj({
        D_("gene", "gene"),
        D_("ncRNA_gene", "gene")
    })),
    
    D_("mirbase", sobj({
        D_("ncRNA_gene", "gene"),
        
        D_("miRNA", "transcript"),
        
        D_("exon", "structure")
    }))
};

void loadGFF(const char *path, const sindex &chridx) {
    GffFile gff;
    AnnotInfo info, sinfo;
    GeneInfo ginfo;
    TranscriptInfo tinfo;
    VariantInfo vinfo;
    FeatureInfo finfo;
    gff.open(path);
    int tblid;
    String tmp;
    stringarray tmps;
    while (gff.next()) {
        auto &data = gff.data();
        if (hs_gff[data.source] && hs_gff[data.source][data.type]) {
            SDictionary row;
            tblid = hs_gff[data.source][data.type];
            if (!chridx.hasKey(data.seqid)) continue;
            row["chromosome"] = chridx[data.seqid];
            row["start"] = data.begin;
            row["end"] = data.end;
            switch(tblid) {
            case GENE:
            {
                row["strand"] = data.strand;
                if (data.attribute.hasKey("gene_id")) gene_id = data.attribute["gene_id"];
                else if (data.attribute.hasKey("ID"))
                    gene_id = data.attribute["ID"].substring(data.attribute["ID"].find(":")+1);
                else {
                    gene_id = "";
                    continue;
                }
                gene_name_idx[gene_id] = tables[GENE].nrow();
                row["gene_ID"] = gene_id;
                
                if (data.attribute.hasKey("Name")) gene_name = data.attribute["Name"];
                else gene_name = gene_id;
                row["NAME"] = gene_name;
                
                int type = 0;
                if (data.attribute.hasKey("biotype")) gene_type = data.attribute["biotype"];
                else gene_type = "";
                if(data.type == "gene") {
                    if (gene_type.match("protein_coding") ||
                        gene_type.match("TEC") ||
                        gene_type.match("TR_") ||
                        gene_type.match("IG_")) type = sbio::PROTEIN_CODING;
                    else if (gene_type.match("lncRNA")) type = LNC_RNA;
                    else if (gene_type.match("ncRNA")) type = NC_RNA;
                    else if (gene_type.match("pseudo")) type = sbio::PSEUDOGENE;
                }
                else if (data.type == "ncRNA_gene") {
                    if (gene_type.match("tRNA")) type = sbio::T_RNA;
                    else if (gene_type.match("rRNA")) type = sbio::R_RNA;
                    else if (gene_type.match("sRNA")) type = sbio::SMALL_RNA;
                    else if (gene_type.match("snRNA")) type = sbio::SN_RNA;
                    else if (gene_type.match("snoRNA")) type = sbio::SNO_RNA;
                    else if (gene_type.match("scRNA")) type = sbio::SC_RNA;
                    else if (gene_type.match("scaRNA")) type = sbio::SCA_RNA;
                    else if (gene_type.match("miRNA")) type = sbio::MI_RNA;
                    else if (gene_type.match("lncRNA")) type = sbio::LNC_RNA;
                    else if (gene_type.match("lincRNA")) type = sbio::LINC_RNA;
                    else if (gene_type.match("antisense")) type = sbio::AS_RNA;
                    else if (gene_type.match("ribozyme")) type = sbio::RIBOZYME;
                    else type = sbio::NON_CODING;
                }
                else if (data.type == "bidirectional_promoter_lncRNA") type = sbio::LNC_RNA;
                else if (data.type == "pseudogene") type = sbio::PSEUDOGENE;
                row["TYPE"] = type;
                if (data.attribute.hasKey("description")) gene_description = data.attribute["description"];
                else gene_description = "";
                row["DESCRIPTION"] = decodeURL(gene_description);
				tables[GENE].addRow(row);
                break;
            }
            else if (tbl_name == "transcript") {
                int gr = -1;
                gene_id = data.attribute["Parent"].substring(data.attribute["Parent"].find(":")+1);
                gr = gene_name_idx[gene_id];
                row["gene_ID"] = gr+1;
                row["NAME"] = data.attribute["Name"];
                
                if (data.type=="mRNA") row["TYPE"] = sbio::M_RNA;
                else if (data.type=="tRNA") row["TYPE"] = sbio::T_RNA;
                else if (data.type=="rRNA") row["TYPE"] = sbio::R_RNA;
                else if (data.type=="ncRNA") row["TYPE"] = sbio::NC_RNA;
                else if (data.type=="lnc_RNA") row["TYPE"] = sbio::LNC_RNA;
                else if (data.type=="scRNA") row["TYPE"] = sbio::SC_RNA;
                else if (data.type=="snRNA") row["TYPE"] = sbio::SN_RNA;
                else if (data.type=="snoRNA") row["TYPE"] = sbio::SNO_RNA;
                else if (data.type=="miRNA") row["TYPE"] = sbio::MI_RNA;
                else if (data.type=="pseudogenic_transcript") row["TYPE"] = sbio::NC_RNA;
                else if (data.type=="transcript") row["TYPE"] = 0;
                else if (data.type=="unconfirmed_transcript") row["TYPE"] = 0;
                else if (data.type=="C_gene_segment") row["TYPE"] = sbio::M_RNA;
                else if (data.type=="D_gene_segment") row["TYPE"] = sbio::M_RNA;
                else if (data.type=="J_gene_segment") row["TYPE"] = sbio::M_RNA;
                else if (data.type=="V_gene_segment") row["TYPE"] = sbio::M_RNA;
                else if (data.type=="three_prime_overlapping_ncrna") row["TYPE"] = sbio::NC_RNA;
                else if (data.type=="vaultRNA_primary_transcript") row["TYPE"] = sbio::VT_RNA;
                
                trans_name_idx[data.attribute["transcript_id"]] = tables[TRANSCRIPT].nrow();
				tables[TRANSCRIPT].addRow(row);
                tables[GENE][gr][trs_col_index_g].add(tables[TRANSCRIPT].addRow());
            }
            else if (tbl_name == "structure") {
                int tr = -1;
                auto trans_id = data.attribute["Parent"].substring(data.attribute["Parent"].find(":")+1);
                tr = trans_name_idx[trans_id];
                
                if (data.type=="CDS") row["TYPE"] = sbio::CDS;
                else if (data.type=="five_prime_UTR") row["TYPE"] = sbio::UTR5;
                else if (data.type=="three_prime_UTR") row["TYPE"] = sbio::UTR3;
                else if (data.type=="exon") row["TYPE"] = sbio::EXON;
                
                row["transcript_ID"] = tr+1;
				tables[STRUCTURE].addRow(row);
                tables[TRANSCRIPT][tr][strct_col_index_t].add(tables[STRUCTURE].nrow());
            }
        }
    }
}
inline void addOrth(const char* path) {



}
extern "C" {
	splugin makeDB(SDataBase& db, const SDictionary& pref, const SeqList& reference) {
        Response res;
        auto nameindex = reference.nameIndex();
        auto suppl = pref["supplementary"] ? pref["supplementary"].parse(";", "=") : sattribute();
        if (suppl["refseq"]) make(pref["orthologue"]);
        loadGFF(pref["source"], nameindex);
        if (suppl["orthologue"]) addOrth(pref["orthologue"]);

        return res;
	}
}
#endif