#include "sapp/scuiapp.h"
#include "sbioinfo/annotation.h"
#include "moirei.h"
using namespace slib;
using namespace slib::smath;
using namespace slib::sio;
using namespace slib::sutil;
using namespace slib::sbio;
using namespace slib::sbio::sio;
using namespace slib::sbio::sutil;

Response moir::gffSummary(moir::Param& param, const SDictionary& pref) {
	GffFile gff;
	GffData* gfd;
	Response res;
	SDictionary dict;
	param.ostream.setStrOStream(res.output);
	gff.open(pref["_args_"][0]);
	while ((gfd = gff.next())) {
		auto& dat = *gfd;
		auto obj = dict[dat.source][dat.type];
		sfor(dat.attribute) {
			if (obj.isArray() && obj.include($_.key())) continue;
			else obj.add($_.key());
		}
	}
	sfor(dict) {
		auto keys = $_.value().keyset();
		keys.sort();
		sforeach(key, keys) {
			param.ostream.print($_.key(), TAB, key, TAB, $_.value()[key]);
		}
	}
	return res;
}

typedef int makeDB(SDataBase& db, const SDictionary& dict, const SeqList& ref);
void initTables(SDataBase& db) {
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
		SColumn("transcript", {D_("type", "integer")}),
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
inline void makeInfoTable(SDataBase &db, const SDictionary& pref) {
	auto table = db["info"];
	table.prepare().insert();
	table.addRecord({ "date", SDate().toString(sstyle::YMD)})
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
	//table.addRow({ "transcript.type", 0x1000, "mRNA" });


	// Genetic region type
	table.addRecord({ "structure.type", sbio::CDS, "CDS", snull })
		.addRecord({ "structure.type", sbio::UTR, "UTR", snull })
		.addRecord({ "structure.type", sbio::UTR5, "5'UTR", snull })
		.addRecord({ "structure.type", sbio::UTR3, "3'UTR", snull })
		.addRecord({ "structure.type", sbio::EXON, "exon", snull })
		.addRecord({ "structure.type", sbio::INTRON, "intron", snull })
		.addRecord({ "structure.type", sbio::SPLICE_SITE, "splice site", snull })
		.addRecord({ "structure.type", sbio::PROCESSING, "processing site (ex. miRNA)", snull });

	// Variant/Mutation type
	table.addRecord({ "variant.type", sbio::SNV, "snv", snull })
		.addRecord({ "variant.type", sbio::MNV, "mnv", snull })
		.addRecord({ "variant.type", sbio::INSERTION, "insertion", snull })
		.addRecord({ "variant.type", sbio::DELETION, "deletion", snull })
		.addRecord({ "variant.type", sbio::DUPLICATION, "duplication", snull })
		.addRecord({ "variant.type", sbio::MULTIPLICATION, "multiplication", snull })
		.addRecord({ "variant.type", sbio::INVERSION, "inversion", snull })
		.addRecord({ "variant.type", sbio::TRANSLOCATION, "translocation", snull });
	
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
	sfori(reference) table.addRecord({ i, reference[i].name, reference[i].attribute["type"], 1, reference[i].length()});
	table.complete();
}
Response moir::makeAnnotDB(moir::Param &par, const SDictionary& pref) {
	try {
		Response res;
		SDataBase db;
		SeqList reference;
		Array<STable> tables;
		reference.open(pref["reference"]);
		if (sfs::exist(pref["output"])) sfs::remove(pref["output"]);
		db.open(pref["output"]);
		initTables(db);
		makeInfoTable(db, pref);
		makeConstTable(db);
		makeChrTable(db, reference);
		sapp::SPlugIn<SDataBase&, const SDictionary&, const SeqList&> plugin(pref["plugin"], "makeDB");
		res = plugin.exec(db, pref, reference);
		return res;
	}
	catch (Exception ex) { ex.print(); return Response(ex); }
}

void getPos(Array<sbpos>& target, String str, sindex& index, AnnotDB& adb){ target.add(RefPos::toPos(str, index)); }
void getPosFromGene(Array<sbpos>& target, String str, sindex& index, AnnotDB& adb) {
	/*
	auto& genes = adb.searchGenes(str);
	if (genes.empty()) target.add(sbpos());
	else target.add(genes[0]);
	*/
}
void getPosFromTranscript(Array<sbpos>& target, String str, sindex& index, AnnotDB& adb) {
	/*
	auto& genes = adb.searchGenes(str);
	if (genes.empty()) target.add(sbpos());
	else {
		sfor(genes) {
			sforeach(transcript, $_.transcripts) {
				target.add(sbpos($_.idx, transcript->begin, transcript->end, $_.dir));
			}
		}
	}
	*/
}
void getPosFromVariant(Array<sbpos>& target, String str, sindex& index, AnnotDB& adb) {
	//auto& vars = adb.searchVariant(str);
	//if (vars.empty()) target.add(sbpos());
	//else target.add(*(vars[0]));
}
void getPosFromMutation(Array<sbpos>& target, String str, sindex& index, AnnotDB& adb) {
	/*
	auto& muts = adb.searchMutations(str);
	if (muts.empty()) target.add(sbpos());
	else target.add(*(muts[0]));
	*/
}
void getPosFromFeature(Array<sbpos>& target, String str, sindex& index, AnnotDB& adb) {
	/*
	auto& ftrs = adb.searchFeature(str);
	if (ftrs.empty()) target.add(sbpos());
	else target.add(*(ftrs[0]));
	*/
}
void printPos(IOStream &stream, sbpos &pos, AnnotDB& adb, const SDictionary& pref) {
	stream << pos.toString(adb.chromosomes) << NL; stream.flush();
}
void printGene(IOStream& stream, sbpos& pos, AnnotDB& adb, const SDictionary &pref) {
	/*
	auto& genes = adb.searchGenes(pos);
	stringarray vals;
	sfor(genes) {
		if (pref["export-geneid"]) vals.add($_.geneid);
		else vals.add($_.name);
	}
	stream << vals << NL; stream.flush();
	*/
}
void printTranscript(IOStream& stream, sbpos& pos, AnnotDB& adb, const SDictionary& pref) {
	/*
	stringarray names;
	auto& genes = adb.searchGenes(pos);
	sfor(genes) {
		sforeach(transcript, $_.transcripts) {
			if (transcript->overlap(pos)) names.add(transcript->name);
		}
	}
	stream << names << NL; stream.flush();
	*/
}
void printVariant(IOStream& stream, sbpos& pos, AnnotDB& adb, const SDictionary& pref) {
	auto& vars = adb.searchVariants(pos);
	stringarray vals;
	sfor(vars) {
		if (pref["exclude-nonsense"]) {
			//if($_->) vals.add($_->gene_id);
		}
		else vals.add($_.name);
	}
	stream << vals << NL; stream.flush();
}
void printMutation(IOStream& stream, sbpos& pos, AnnotDB& adb, const SDictionary& pref) {
	auto& muts = adb.searchMutations(pos);
	stringarray vals;
	sfor(muts) {
		if (pref["exclude-nonsense"]) {
			//if($_->) vals.add($_->gene_id);
		}
		else vals.add($_.name);
	}
	stream << vals << NL; stream.flush();
}
void printFeature(IOStream& stream, sbpos& pos, AnnotDB& adb, const SDictionary& pref) {
	//auto& ftrs = adb.searchVariants(pos);
	stringarray vals;
	//sfor(ftrs) { vals.add($_->name); }
	stream << vals << NL; stream.flush();
}
Response moir::bioAnnot(moir::Param& param, const SDictionary& pref) {
	Response res;
	Array<sbpos> target;
	AnnotDB adb(pref["annotdb"]);
	IOStream outstream;
	String ln;
	SFile of;
	if (pref["output"]) {
		of.open(pref["output"], MAKE);
		outstream.setFileOStream(of);
	}
	else outstream.setStrOStream(res.output);
	if (pref["from"] == pref["to"]) {
		if (pref["input"]) {
			SFile in(pref["input"]);
			while (in) {
				in.readLine(ln);
				outstream << ln << NL;
			}
		}
		else outstream << pref["_args_"];
		outstream.flush();
	}
	else {
		SFunction<void, Array<sbpos>&, String, sindex&, AnnotDB&> importer;
		SFunction<void, IOStream&, sbpos&, AnnotDB&, const SDictionary&> exporter;
		if (pref["from"] == "pos") importer = getPos;
		else if (pref["from"] == "gene") {
			importer = getPosFromGene;
			adb.load({ "gene" });
		}
		else if (pref["from"] == "transcript") importer = getPosFromTranscript;
		else if (pref["from"] == "variant") importer = getPosFromVariant;
		else if (pref["from"] == "mutation") importer = getPosFromMutation;
		else if (pref["from"] == "feature") importer = getPosFromFeature;
		if (pref["to"] == "pos") exporter = printPos;
		else if (pref["to"] == "gene") exporter = printGene;
		else if (pref["to"] == "transcript") exporter = printTranscript;
		else if (pref["to"] == "variant") exporter = printVariant;
		else if (pref["to"] == "mutation") exporter = printMutation;
		else if (pref["to"] == "feature") exporter = printFeature;
		try {
			if (pref["input"]) {
				String ln;
				SFile f(pref["input"]);
				while (f) {
					f.readLine(ln);
					importer(target, ln, adb.chrindex, adb);
				}
			}
			else {
				sfor(pref["_args_"]) importer(target, $_.toString(), adb.chrindex, adb);
			}
			sfor(target) exporter(outstream, $_, adb, pref);
		}
		catch (Exception ex) {
			ex.print();
		}
	}
	return res;
}


