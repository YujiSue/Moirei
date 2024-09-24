#include "moirei.h"
Response& Moirei::convertGenome(const SDictionary& pref) {
	try {
		param.setPref(pref);
		String refpath, source = pref.hasKey("source") ? (const char*)pref["source"] : "";
		// Download if input is FTP or HTTP(S) URL
		if (pref["input"].match(REG("/^https?:|ftp:/"))) {
			if (pref["input"].match("ensembl") && source.empty()) source = "ensembl";
			else if (pref["input"].match("refseq") && source.empty()) source = "refseq";
			auto res = download(pref["input"], { D_("outdir", param.outdir) });
			if (res.code) throw Exception(res);
			refpath = res.output;
		}
		// Use local file
		else refpath = pref["input"];
		/*
		// Check file is archived
		auto ext = sstr::toLower(sfs::extension(refpath));
		if (ext == "zip" || ext == "gz") {}
		*/
		// Make index
		Fasta fa;
		SeqList reference;
		auto seqtype = DNA_SEQ4 | sseq::MASKED;
		SPrint("Open reference file '", refpath, "'.");
		fa.open(refpath, seqtype);
		SPrint("Making index...");
		fa.makeIndex();
		sveci available;
		available.resize(fa.count(), 1);
		// Set conversion target
		int target = 0;
		if (pref["lg"] == "ALL") target = 0xFF;
		else {
			if (pref["lg"].match("CHR")) target |= (int)LG_TYPE::CHROMOSOME;
			if (pref["lg"].match("M")) target |= (int)LG_TYPE::MT_GENOME;
			if (pref["lg"].match("PLT")) target |= (int)LG_TYPE::PL_GENOME;
			if (pref["lg"].match("PLM")) target |= (int)LG_TYPE::PLASMID;
		}
		// (Optional) Load sequence attribute,
		SDictionary attr;
		if (pref.hasKey("attribute") && sfs::exist(pref["attribute"])) {
			SPrint("Load attribute.");
			SFile f(pref["attribute"]);
			if (source == "refseq") {
				while (f) {
					f.readLine(param.ln);
					if (param.ln.empty() || param.ln[0] == '#') continue;
					auto vals = param.ln.split(TAB);
					attr[vals[6]] = {
						D_("ori", "genbank"),
						D_("name", vals[0]),
						D_("type", (int)LG_TYPE::MISC_LG)
					};
					if (vals[1] == "assembled-molecule") {
						if (vals[3] == "Chromosome") attr[vals[6]]["type"] = (int)LG_TYPE::CHROMOSOME;
						else if (vals[3] == "Mitochondrion") attr[vals[6]]["type"] = (int)LG_TYPE::MT_GENOME;
					}
					//attr[vals[6]]["form"] = 0;
				}
			}
			else {
				while (f) {
					f.readLine(param.ln);
					if (param.ln.empty() || param.ln[0] == '#') continue;
					auto vals = param.ln.split(TAB);
					attr[vals[0]] = {
						D_("ori", source),
						D_("name", vals[1]),
						D_("type", (int)LG_TYPE::MISC_LG)
					};
					vals[2] = sstr::toLower(vals[2]);
					if (vals[2].match("auto")) attr[vals[0]]["type"] = (int)LG_TYPE::AUTOSOME;
					else if (vals[2].match("sex")) attr[vals[0]]["type"] = (int)LG_TYPE::SEX_CHROMOSOME;
					else if (vals[2].match("some")) attr[vals[0]]["type"] = (int)LG_TYPE::CHROMOSOME;
					else if (vals[2].match("mt") || vals[2].match("mito")) attr[vals[0]]["type"] = (int)LG_TYPE::MT_GENOME;
					else if (vals[2].match("plastid")) attr[vals[0]]["type"] = (int)LG_TYPE::PL_GENOME;
					else if (vals[2].match("plasmid")) attr[vals[0]]["type"] = (int)LG_TYPE::PLASMID;
					//attr[vals[0]]["form"] = vals[3] == "c" ? (int)sbio::sseq::CIRCULAR : 0;
				}
			}
			sfor2(fa.titles, available) {
				if (attr.hasKey($_1) && (attr[$_1]["type"].ubyteValue() & target)) $_2 = 1;
				else $_2 = 0;
			}
		}
		// Set general info.
		if (source.size()) reference.attribute["source"] = source;
		if (pref.hasKey("species")) reference.attribute["species"] = pref["species"];
		if (pref.hasKey("ver")) reference.attribute["version"] = pref["ver"];
		// 
		reference.resize(sstat::sum(available));
		if (!param.output.endWith(".bin")) param.output += ".bin";
		param.response.output = param.output;
		// Conversion
		SPrint("Started conversion.");
		auto it = reference.begin();
		sforin(i, 0, fa.count()) {
			if (!available[i]) continue;
			fa.setIndex(i);
			fa >> $_;
			if (attr[$_.name] && attr[$_.name]["type"]) $_.attribute["type"] = attr[$_.name]["type"];
			else $_.attribute["type"] = (int)LG_TYPE::MISC_LG;
			if (pref["rename"]) {
				$_.attribute.set(attr[$_.name]["ori"], $_.name);
				$_.name = attr[$_.name]["name"];
			}
			$NEXT;
		}
		SPrint("Completed.");
		// Save to
		SPrint("Save to '", param.output, "'.");
		reference.save(param.output);
	}
	catch (Exception ex) {
		ex.print();
		param.response = ex;
	}
	return param.response;
}
