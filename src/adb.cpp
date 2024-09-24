#include "moirei.h"
Response& Moirei::gffSummary(const SDictionary& pref) {
	GffFile gff;
	GffData* gfd;
	SDictionary dict;
	size_t size = 0, current = 0;
	bool load = true;
	if (pref["verbose"]) SPrint("Open GFF file '", pref["_args_"][0], "'.");
	try {
		gff.open(pref["_args_"][0]);
		size = gff.size();
		while ((gfd = gff.next())) {
			auto& dat = *gfd;
			auto obj = dict[dat.source][dat.type];
			sfor(dat.attribute) {
				if (obj.isArray() && obj.include($_.key())) continue;
				else obj.add($_.key());
			}
		}
		if (param.output.empty()) param.ostream.setStrOStream(param.response.output);
		else {
			param.ofile.open(param.output, MAKE);
			param.ostream.setFileOStream(param.ofile);
			param.response.output = param.output;
		}
		param.ostream.print(S("=") * 60);
		param.ostream.print("Source");
		param.ostream.print(SP * 3, "|-  Type : [ Attribute keys... ]");
		param.ostream.print(S("=") * 60);
		auto keys1 = dict.keyset();
		keys1.sort();
		sforeach(key1, keys1) {
			param.ostream.print(S("- ") * 30);
			param.ostream.print(key1);
			auto keys2 = dict[key1].keyset();
			keys2.sort();
			sforeach(key2, keys2) {
				param.ostream.print(SP * 3, "|-  ", key2, " : ", dict[key1][key2]);
			}
		}
		param.ostream.print(S("- ") * 30);
		if (param.response.output) SPrint(param.response.output);
	}
	catch (Exception ex) {
		ex.print();
		param.response = ex;
	}
	return param.response;
}
typedef int makeDB(SDataBase& db, const SDictionary& dict, const SeqList& ref);
inline void initTables(SDataBase& db) {
	// info
	db.create("info", {
		SColumn("key", {D_("type", "string")}),
		SColumn("value", {D_("type", "string")})
		});
	// const
	db.create("const", {
		SColumn("target", {D_("type", "string")}),
		SColumn("value", {D_("type", "integer")}),
		SColumn("label", {D_("type", "string")}),
		SColumn("description", {D_("type", "string")})
		});
	// chromosome
	db.create("chromosome", {
		SColumn("id", {D_("type", "integer"), D_("rule", "key")}),
		SColumn("name", {D_("type", "string")}),
		SColumn("type", {D_("type", "integer")}),
		SColumn("start", {D_("type", "integer")}),
		SColumn("end", {D_("type", "integer")})
		});
	db.create("contig", {
		SColumn("name", {D_("type", "string")}),
		SColumn("chromosome", {D_("type", "integer")}),
		SColumn("start", {D_("type", "integer")}),
		SColumn("end", {D_("type", "integer")})
		});
	db.create("gene", {
		SColumn("id", {D_("type", "integer"), D_("rule", "key")}),
		SColumn("name", {D_("type", "string")}),
		SColumn("type", {D_("type", "integer")}),
		SColumn("chromosome", {D_("type", "integer")}),
		SColumn("start", {D_("type", "integer")}),
		SColumn("end", {D_("type", "integer")}),
		SColumn("strand", {D_("type", "bool")}),
		SColumn("gid", {D_("type", "string")}),
		SColumn("description", {D_("type", "string")}),
		SColumn("attribute", {D_("type", "dict"), D_("format", "json")})
		});
	db.create("synonym", {
		SColumn("gene", {D_("type", "integer")}),
		SColumn("name", {D_("type", "string")})
		});
	db.create("transcript", {
		SColumn("id", {D_("type", "integer"), D_("rule", "key")}),
		SColumn("name", {D_("type", "string")}),
		SColumn("type", {D_("type", "integer")}),
		SColumn("gene", {D_("type", "integer")}),
		SColumn("start", {D_("type", "integer")}),
		SColumn("end", {D_("type", "integer")})
		});
	db.create("structure", {
		SColumn("type", {D_("type", "integer")}),
		SColumn("transcript", {D_("type", "integer")}),
		SColumn("start", {D_("type", "integer")}),
		SColumn("end", {D_("type", "integer")})
		});
	db.create("mutation", {
		SColumn("id", {D_("type", "integer"), D_("rule", "key")}),
		SColumn("name", {D_("type", "string")}),
		SColumn("type", {D_("type", "integer")}),
		SColumn("chromosome", {D_("type", "integer")}),
		SColumn("start", {D_("type", "integer")}),
		SColumn("end", {D_("type", "integer")}),
		SColumn("vid", {D_("type", "string")}),
		SColumn("attribute", {D_("type", "dict"), D_("format", "json")})
		});
	db.create("variant", {
		SColumn("id", {D_("type", "integer"), D_("rule", "key")}),
		SColumn("name", {D_("type", "string")}),
		SColumn("type", {D_("type", "integer")}),
		SColumn("chromosome", {D_("type", "integer")}),
		SColumn("start", {D_("type", "integer")}),
		SColumn("end", {D_("type", "integer")}),
		SColumn("vid", {D_("type", "string")}),
		SColumn("attribute", {D_("type", "dict"), D_("format", "json")})
		});
	db.create("feature", {
		SColumn("name", {D_("type", "string")}),
		SColumn("type", {D_("type", "integer")}),
		SColumn("chromosome", {D_("type", "integer")}),
		SColumn("start", {D_("type", "integer")}),
		SColumn("end", {D_("type", "integer")}),
		SColumn("strand", {D_("type", "bool")})
		});
	db.create("protein", {
		SColumn("id", {D_("type", "integer"), D_("rule", "key")}),
		SColumn("name", {D_("type", "string")}),
		SColumn("gene", {D_("type", "integer")}),
		SColumn("start", {D_("type", "integer")}),
		SColumn("end", {D_("type", "integer")})
		});
	db.create("motif", {
		SColumn("name", {D_("type", "string")}),
		SColumn("type", {D_("type", "integer")}),
		SColumn("protein", {D_("type", "integer")}),
		SColumn("start", {D_("type", "integer")}),
		SColumn("end", {D_("type", "integer")}),
		SColumn("mid", {D_("type", "string")}),
		SColumn("source", {D_("type", "string")})
		});
	db.create("xref", {
		SColumn("type", {D_("type", "integer")}),
		SColumn("dbid", {D_("type", "integer")}),
		SColumn("ref", {D_("type", "string")}),
		SColumn("refid", {D_("type", "string")})
		});
}
inline void makeInfoTable(SDataBase& db, const SDictionary& pref) {
	auto table = db["info"];
	table.prepare().insert();
	table.addRecord({ "date", SDate().toString(sstyle::YMD) })
		.addRecord({ "species", pref["species"] })
		.addRecord({ "source", pref["source"] })
		.addRecord({ "version", pref["dbver"] });
	table.complete();
}
inline void makeConstTable(SDataBase& db) {
	auto table = db["const"];
	table.prepare().insert();
	// Linkage group type
	table.addRecord({ "chromosome.type", (int)sbio::LG_TYPE::CHROMOSOME, "chromosome", snull })
		.addRecord({ "chromosome.type", (int)sbio::LG_TYPE::AUTOSOME, "autosome", snull })
		.addRecord({ "chromosome.type", (int)sbio::LG_TYPE::SEX_CHROMOSOME, "sex chromosome", snull })
		.addRecord({ "chromosome.type", (int)sbio::LG_TYPE::MT_GENOME, "mitochondria genome", snull })
		.addRecord({ "chromosome.type", (int)sbio::LG_TYPE::PL_GENOME, "plastid genome", snull })
		.addRecord({ "chromosome.type", (int)sbio::LG_TYPE::PLASMID, "plasmid", snull })
		.addRecord({ "chromosome.type", (int)sbio::LG_TYPE::MISC_LG, "misc linkage group", snull });
	// Gene type
	table.addRecord({ "gene.type", (int)sbio::GENE_TYPE::PROTEIN_CODING,"protein coding", snull })
		.addRecord({ "gene.type", (int)sbio::GENE_TYPE::NON_CODING, "non-coding", snull })
		.addRecord({ "gene.type", (int)sbio::GENE_TYPE::PSEUDO_GENE, "pseudogene", snull })
		.addRecord({ "gene.type", (int)sbio::GENE_TYPE::TRANSPOSON, "transposon",  snull })
		.addRecord({ "gene.type", (int)sbio::GENE_TYPE::MISC_GENE, "gene",  snull })
		.addRecord({ "gene.type", (int)sbio::GENE_TYPE::UNAVAILABLE, "dead",  snull });

	// Transcript type
	table.addRecord({ "transcript.type", (int)sbio::TRANSCRIPT_TYPE::M_RNA , "mRNA", snull })
		.addRecord({ "transcript.type", (int)sbio::TRANSCRIPT_TYPE::PSEUDO_GENE_TRANSCRIPT , "pseudogene transcript", snull })
		.addRecord({ "transcript.type", (int)sbio::TRANSCRIPT_TYPE::T_RNA , "tRNA", snull })
		.addRecord({ "transcript.type", (int)sbio::TRANSCRIPT_TYPE::R_RNA , "rRNA", snull })
		.addRecord({ "transcript.type", (int)sbio::TRANSCRIPT_TYPE::NC_RNA , "ncRNA", snull })
		.addRecord({ "transcript.type", (int)sbio::TRANSCRIPT_TYPE::LNC_RNA , "long ncRNA", snull })
		.addRecord({ "transcript.type", (int)sbio::TRANSCRIPT_TYPE::LINC_RNA , "long intergenic ncRNA", snull })
		.addRecord({ "transcript.type", (int)sbio::TRANSCRIPT_TYPE::SC_RNA , "small conditional RNA", snull })
		.addRecord({ "transcript.type", (int)sbio::TRANSCRIPT_TYPE::SN_RNA , "small nuclear RNA", snull })
		.addRecord({ "transcript.type", (int)sbio::TRANSCRIPT_TYPE::SNO_RNA , "small nucleolar RNA", snull })
		.addRecord({ "transcript.type", (int)sbio::TRANSCRIPT_TYPE::TEROMERASE_RNA , "telomerase RNA" })
		.addRecord({ "transcript.type", (int)sbio::TRANSCRIPT_TYPE::Y_RNA , "Y RNA", snull })
		.addRecord({ "transcript.type", (int)sbio::TRANSCRIPT_TYPE::VT_RNA , "vault RNA" })
		.addRecord({ "transcript.type", (int)sbio::TRANSCRIPT_TYPE::RNASE_P , "RNaseP" })
		.addRecord({ "transcript.type", (int)sbio::TRANSCRIPT_TYPE::MI_RNA , "miRNA" })
		.addRecord({ "transcript.type", (int)sbio::TRANSCRIPT_TYPE::PI_RNA , "piRNA" })
		.addRecord({ "transcript.type", (int)sbio::TRANSCRIPT_TYPE::AS_RNA , "antisense RNA" });

	// Genetic region type
	table.addRecord({ "structure.type", sbio::CDS, "CDS", snull })
		.addRecord({ "structure.type", sbio::UTR, "UTR", snull })
		.addRecord({ "structure.type", sbio::UTR5, "5'UTR", snull })
		.addRecord({ "structure.type", sbio::UTR3, "3'UTR", snull })
		.addRecord({ "structure.type", sbio::EXON, "exon", snull })
		.addRecord({ "structure.type", sbio::INTRON, "intron", snull })
		.addRecord({ "structure.type", sbio::SPLICE_SITE, "splice site", snull })
		.addRecord({ "structure.type", sbio::PROCESSING, "processing site (ex. miRNA)", snull });

	// Variant,Mutation type
	table.addRecord({ "variant.type", sbio::SNV, "snv", snull })
		.addRecord({ "variant.type", sbio::MNV, "mnv", snull })
		.addRecord({ "variant.type", sbio::INSERTION, "insertion", snull })
		.addRecord({ "variant.type", sbio::DELETION, "deletion", snull })
		.addRecord({ "variant.type", sbio::DUPLICATION, "duplication", snull })
		.addRecord({ "variant.type", sbio::MULTIPLICATION, "multiplication", snull })
		.addRecord({ "variant.type", sbio::INVERSION, "inversion", snull })
		.addRecord({ "variant.type", sbio::TRANSLOCATION, "translocation", snull })
		.addRecord({ "variant.type", (sbio::CDS << 8), "CDS variant", snull })
		.addRecord({ "variant.type", (sbio::EXON << 8), "exon variant", snull })
		.addRecord({ "variant.type", (sbio::UTR << 8), "UTR variant", snull })
		.addRecord({ "variant.type", (sbio::INTRON << 8), "intron variant", snull })
		.addRecord({ "variant.type", (sbio::SPLICE_SITE << 8), "splice site variant", snull })
		.addRecord({ "variant.type", (sbio::PROCESSING << 8), "RNA processing site variant", snull })
		.addRecord({ "variant.type", (sbio::MISSENSE << 16), "missense variant", snull })
		.addRecord({ "variant.type", (sbio::NONSENSE << 16), "nonsense variant", snull })
		.addRecord({ "variant.type", (sbio::SUBSTITUTION << 16), "base substitution", snull })
		.addRecord({ "variant.type", (sbio::FRAME_SHIFT << 16), "frame-shift variant", snull })
		.addRecord({ "variant.type", (sbio::IN_FRAME << 16), "in-frame variant", snull })
		.addRecord({ "variant.type", (sbio::INDEL << 16), "short indel", snull });

	// Feature flag
	//table.addRow({ snull, "feature.type", 0x00, "" });

	// Cross reference
	table.addRecord({ "xref.type", 1, "gene", snull })
		.addRecord({ "xref.type", 2, "transcript", snull })
		.addRecord({ "xref.type", 3, "protein", snull });
	table.complete();
}
inline void makeChrTable(SDataBase& db, const SeqList& reference) {
	auto table = db["chromosome"];
	table.prepare().insert();
	sfori(reference) table.addRecord({ i, reference[i].name, reference[i].attribute["type"], 1, reference[i].length() });
	table.complete();
}
Response& Moirei::makeAnnotDB(const SDictionary& pref) {
	try {
		param.setPref(pref);
		//
		SDataBase db;
		SeqList reference;
		Array<STable> tables;
		reference.open(pref["reference"]);
		if (sfs::exist(param.output)) sfs::remove(pref["output"]);
		if (!param.output.endWith(".db")) param.output += ".db";
		param.response.output = param.output;
		//
		db.open(param.output);
		//
		initTables(db);
		makeInfoTable(db, pref);
		makeConstTable(db);
		makeChrTable(db, reference);
		sapp::SPlugIn<SDataBase&, const SDictionary&, const SeqList&> plugin(pref["plugin"], "makeDB");
		param.response.code = plugin.exec(db, pref, reference);
	}
	catch (Exception ex) { 
		ex.print(); 
		return param.response = ex;
	}
	return param.response;
}