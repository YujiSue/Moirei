#include "moirei.h"
using namespace slib;
using namespace slib::sio;
using namespace slib::sutil;
using namespace slib::sbio;
using namespace moir;

Response moir::convertGenome(moir::Param& param, const SDictionary& pref) {
	try {
		String refpath, source = pref["source"] ? (const char *)pref["source"] : "";
		Response res;
		// Download if input is FTP or HTTP(S) URL
		if (pref["input"].match(REG_(/^https?:|ftp:/))) {
			if (pref["input"].match("ensembl") && source.empty()) source = "ensembl";
			else if (pref["input"].match("refseq") && source.empty()) source = "refseq";
			res = download(pref["input"], { D_("outdir", pref["outdir"]) });
			if (res.code) throw Exception(res);
			refpath = res.output;
		}
		// Use local file
		else refpath = pref["input"];
		// Check file is archived 
		auto ext = sstr::toLower(sfs::extension(refpath));
		if (ext == "zip" || ext == "gz") {}
		// Make index of fasta
		Fasta fa;
		SeqList reference;
		auto seqtype = DNA_SEQ4 | sseq::MASKED;
		fa.open(refpath, seqtype);
		fa.makeIndex();
		sveci available;
		available.resize(fa.count(), 1);
		//
		int target = 0;
		if (pref["target-lg"] == "ALL") target = 0xFF;
		else {
			if (pref["target-lg"].match("CHR")) target |= (int)LG_TYPE::CHROMOSOME;
			if (pref["target-lg"].match("M")) target |= (int)LG_TYPE::MT_GENOME;
			if (pref["target-lg"].match("PLT")) target |= (int)LG_TYPE::PL_GENOME;
			if (pref["target-lg"].match("PLM")) target |= (int)LG_TYPE::PLASMID;
		}
		// (Optional) Sequence attribute to change names or set flag of linkage group type,
		SDictionary attr;		
		if (pref["attribute-file"] && sfs::exist(pref["attribute-file"])) {
			SFile f(pref["attribute-file"]);
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
		// Set attributes
		if (source.size()) reference.attribute["source"] = source;
		if (pref["species"]) reference.attribute["species"] = pref["species"];
		if (pref["refver"]) reference.attribute["version"] = pref["refver"];
		// 
		reference.resize(sstat::sum(available));
		res.output = sfs::joinPath(pref["outdir"], pref["file-name"] + ".bin");
		// Conversion
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
		// Save to
		reference.save(res.output);
		return res;
	}
	catch (Exception ex) {
		ex.print();
		return Response(ex);
	}
}
Response moir::extractGenome(moir::Param& param, const SDictionary& pref) {
	try {
		Response res;
		RefPos pos;
		SeqList reference, seqs;
		reference.load(pref["reference"]);
		auto nameindex = reference.nameIndex();
		bool annot = false;
		AnnotDB db;
		if (pref.hasKey("output")) {
			param.ofile.open(pref["output"], MAKE);
			if (param.oformat == "auto") param.oformat = sfs::extension(pref["output"]);
			param.ostream.setFileOStream(param.ofile);
		}
		else {
			if (param.oformat == "auto") param.oformat = "fa";
			param.ostream.setStrOStream(res.output);
		}
		if (pref.hasKey("annotation")) annot = true;
		if (annot) db.open(pref["annotdb"]);
		auto sites = pref["_args_"];
		seqs.resize(sites.size());
		sfor2(sites, seqs) {
			pos = RefPos::toPos($_1, nameindex);
			$_2.setSeqAs(reference.raw(pos), sseq::DNA);
			$_2.name = $_1;
			if (annot) {
				auto &genes = db.getGenes(pos);
			}
		}
		sbio::sio::writeSeqs(param.ostream, seqs, param.oformat);
		return res;
	}
	catch (Exception ex) {
		ex.print();
		return Response(ex);
	}
}
Response moir::referenceSummary(moir::Param& param, const SDictionary& pref) {
	try {
		Response res;
		SeqList reference;
		reference.open(pref["reference"]);
		AnnotDB adb;
		size_t count[3] = { 0, 0, 0 };
		adb.open(pref["annotdb"]);
		auto geneTbl = adb["gene"];
		count[0] = geneTbl.where("type&" + S((int)GENE_TYPE::PROTEIN_CODING)).count(); geneTbl.reset();
		count[1] = geneTbl.where("type&" + S((int)GENE_TYPE::NON_CODING)).count(); geneTbl.reset();
		count[2] = geneTbl.where("type!=0").count();
		if (pref.hasKey("output")) {
			param.ofile.open(pref["output"], MAKE);
			param.ostream.setFileOStream(param.ofile);
		}
		else param.ostream.setStrOStream(res.output);
		param.ostream.print(S_(= ) * 60);
		param.ostream.print(sstr::fill("Count:", ' ', 10, DIRECTION::TAIL), reference.size());
		param.ostream.print(sstr::fill("Source:", ' ', 10, DIRECTION::TAIL), reference.attribute["source"]);
		param.ostream.print(sstr::fill("Species:", ' ', 10, DIRECTION::TAIL), reference.attribute.hasKey("species") ? reference.attribute["species"] : "Unknown");
		param.ostream.print(sstr::fill("Version:", ' ', 10, DIRECTION::TAIL), reference.attribute.hasKey("version") ? reference.attribute["version"] : "Unknown");
		param.ostream.print(NL);
		//
		param.ostream.print(S_(-) * 40);
		param.ostream.print("|", sstr::fill("Name", ' ', 12, DIRECTION::BI),
			"|", sstr::fill("Size (bp)", ' ', 12, DIRECTION::BI),
			"|", sstr::fill("Mask (bp)", ' ', 12, DIRECTION::BI),
			"|");
		param.ostream.print(S_(-) * 40);
		sfor(reference) {
			param.ostream.print("|", sstr::fill($_.name, ' ', 12, DIRECTION::TAIL),
				"|", sstr::fill(String($_.length()), ' ', 12, DIRECTION::HEAD),
				"|", sstr::fill(String($_.mask.length(true)), ' ', 12, DIRECTION::HEAD),
				"|");
		}
		param.ostream.print(S_(-) * 40);
		param.ostream.print(NL);
		param.ostream.print(S_(-) * 40);
		//
		param.ostream.print("|", sstr::fill("Annotated", ' ', 25, DIRECTION::BI),
			"|", sstr::fill("Count", ' ', 12, DIRECTION::BI),
			"|");
		param.ostream.print(S_(-) * 40);
		param.ostream.print("|", sstr::fill("Genes", ' ', 25, DIRECTION::TAIL),
			"|", sstr::fill(String(count[2]), ' ', 12, DIRECTION::BI),
			"|");
		param.ostream.print("|", SP * 5, sstr::fill("Protein coding", ' ', 20, DIRECTION::TAIL),
			"|", sstr::fill(String(count[0]), ' ', 12, DIRECTION::BI),
			"|");
		param.ostream.print("|", SP * 5, sstr::fill("Non-coding", ' ', 20, DIRECTION::TAIL),
			"|", sstr::fill(String(count[1]), ' ', 12, DIRECTION::BI),
			"|");
		param.ostream.print("|", SP * 5, sstr::fill("Others", ' ', 20, DIRECTION::TAIL),
			"|", sstr::fill(String(count[2]-count[0]-count[1]), ' ', 12, DIRECTION::BI),
			"|");
		param.ostream.print(S_(-) * 40);
		param.ostream.print(S_(= ) * 60);
		return res;
	}
	catch (Exception ex) {
		ex.print();
		return Response(ex);
	}
}
Response moir::countGCRatio(moir::Param& param, const SDictionary& pref) {
	try {
		Response res;
		SeqList reference;
		ubytearray seq;
		reference.load(pref["reference"]);
		if (pref.hasKey("output")) param.ofile.open(pref["output"], MAKE);
		else param.ofile.open(sfs::joinPath(sfs::splitPath(pref["reference"]).first, sfs::fileName(pref["reference"], false) + "_gc.bin"), MAKE);
		//
		param.ofile.writeInt((int)reference.size());
		int bin = pref["bin"], n = bin / 4;
		bin = n * 4;
		seq.resize(bin);
		//
		param.ofile.writeInt(bin);
		sfori(reference) {
			//
			param.ofile.writeInt(i);
			//
			param.ofile.writeInt((int)((reference[i].size() - 1) / n + 1));
			auto bp = reference[i].data();
			size_t current = 0, end = reference[i].size() - n;
			while (current < end) {
				sdna::expand4(bp, 0, bin, seq.data());
				//
				param.ofile.writeInt((int)sna::gcCount(seq));
				bp += n; current += n;
			}
			seq.resize(reference[i].size() - current);
			sdna::expand4(bp, 0, seq.size(), seq.data());
			//
			param.ofile.writeInt((int)sna::gcCount(seq));
		}
		return res;
	}
	catch (Exception ex) {
		ex.print();
		return Response(ex);
	}

}