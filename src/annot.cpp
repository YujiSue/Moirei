#include "moirei.h"

void getPos(Array<sbpos>& target, String str, AnnotDB& adb, const SDictionary &pref) { 
	target.add(RefPos::toPos(str, adb.chrindex)); 
}

void getPosFromGene(Array<sbpos>& target, String str, AnnotDB& adb, const SDictionary& pref) {
	sobj seachopt = { D_("search-synonym", true), D_("search-xref", true) }, 
		getopt = { D_("select", "id,chromosome,start,end,strand,gid,name") };
	auto& gidx = adb.searchGenes(str, MATCH::EXACT, seachopt);
	if (gidx.size()) target.add(adb.geneInfo(gidx[0], getopt));
}
void getPosFromTranscript(Array<sbpos>& target, String str, AnnotDB& adb, const SDictionary& pref) {
	sobj opt = { D_("select", "gene.chromosome,transcript.start,transcript.end") };
	auto rna = adb.getTranscripts(str, MATCH::EXACT, opt);
	if (rna.size()) target.add(rna[0]);
	else {
		auto rnas = adb.getTranscripts(str, MATCH::PARTIAL, opt);
		sfor(rnas) target.add($_);
	}
}
void getPosFromVariant(Array<sbpos>& target, String str, AnnotDB& adb, const SDictionary& pref) {
	sobj getopt = { D_("select", "chromosome,start,end,vid,name") };
	auto& muts = adb.getMutations(str, MATCH::EXACT, getopt);
	if (muts.size()) {
		target.add(muts[0]); return;
	}
	auto &vars = adb.getVariants(str, MATCH::EXACT, getopt);
	if (vars.size()) {
		target.add(vars[0]); return;
	}
}

void setPos(SArray& array, sbpos& pos, AnnotDB& adb, const SDictionary& pref) {
	array.add(pos.toString(adb.chromosomes));
}
void setGene(SArray& array, sbpos& pos, AnnotDB& adb, const SDictionary& pref) {
	sobj opt = { D_("search-synonym", true), D_("search-xref", true), D_("search-description", true) };
	if (pref.hasKey("export-geneid") && pref["export-geneid"]) opt["select"] = "gid";
	else opt["select"] = "name";
	auto& genes = adb.getGenes(pos, opt);
	sfor(genes) array.add($_.name);
}
void setTranscript(SArray& array, sbpos& pos, AnnotDB& adb, const SDictionary& pref) {
	sobj opt = { D_("select", "name") };
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
	auto& transcripts = adb.getTranscripts(pos, opt);
	sfor(transcripts) array.add($_.name);
}
void setVariant(SArray& array, sbpos& pos, AnnotDB& adb, const SDictionary& pref) {
	sobj opt = { D_("select", "name") };
	if (pref.hasKey("mutation-type")) opt["type"] = pref["mutation-type"];
	if (pref.hasKey("cds-only")) opt["cds-only"] = pref["cds-only"];
	auto& muts = adb.getMutations(pos, opt);
	sfor(muts) { array.add($_.name); }
	auto& vars = adb.getVariants(pos, opt);
	sfor(vars) { array.add($_.name); }
}
void setFeature(SArray& array, sbpos& pos, AnnotDB& adb, const SDictionary& pref) {
	sobj opt = { D_("select", "name") };
	if (pref.hasKey("feature-type")) opt["type"] = pref["feature-type"];
	auto& ftrs = adb.getFeatures(pos, opt);
	sfor(ftrs) { array.add($_.name); }
}
Response& Moirei::bioAnnot(const SDictionary& pref) {
	try {
		param.setPref(pref);
		//
		AnnotDB adb(pref["annotdb"]);
		stringarray inputs;
		SArray result;
		//
		if (pref.hasKey("input")) {
			SFile f(pref["input"]);
			while (f) {
				f.readLine(param.ln);
				inputs.add(sstr::trim(param.ln));
			}
		}
		else {
			sfor(pref["_args_"]) inputs.add($_);
		}
		//
		if (pref["from"] == pref["to"]) {
			/*
			*/
			sfor(inputs) {
				/*
				*/
				result.add(sobj({ D_("query", $_), D_("out", SArray({$_})) }));
			}
		}
		else {
			Array<sbpos> target;
			//
			SFunction<void, Array<sbpos>&, String, AnnotDB&, const SDictionary&> importer;
			SFunction<void, SArray&, sbpos&, AnnotDB&, const SDictionary&> exporter;
			if (pref["from"] == "pos") importer = getPos;
			else if (pref["from"] == "gene") importer = getPosFromGene;
			else if (pref["from"] == "transcript") importer = getPosFromTranscript;
			else if (pref["from"] == "variant" || pref["from"] == "mutation") importer = getPosFromVariant;
			if (pref["to"] == "pos") exporter = setPos;
			else if (pref["to"] == "gene") exporter = setGene;
			else if (pref["to"] == "transcript") exporter = setTranscript;
			else if (pref["to"] == "variant") exporter = setVariant;
			else if (pref["to"] == "feature") exporter = setFeature;
			//
			sfor(inputs) {
				importer(target, $_, adb, pref);
				result.add({ D_("query", $_), D_("out", SArray()) });
				sforeach(pos, target) {
					exporter(result[-1]["out"].array(), pos, adb, pref);
				}
			}
		}
		if (param.oformat == "_obj_") param.response.attribute["result"] = result;
		else {
			if (param.output.size()) {
				param.ofile.open(param.output, MAKE);
				param.ostream.setFileOStream(param.ofile);
			}
			else param.ostream.setStrOStream(param.response.output);
			param.ostream << "Query" << TAB << "Result" << LF;
			sfor(result) {
				param.ostream << $_["query"] << TAB << ($_["out"].empty() ? "Not found." : $_["out"].toString("csv")) << LF;
			}
			if (param.response.output.size() && !pref["silent"]) SPrint(param.response.output);
		}
	}
	catch (Exception ex) {
		ex.print();
		param.response = ex;
	}
	return param.response;
}
