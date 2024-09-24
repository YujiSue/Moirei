#include "moirei.h"
Response& Moirei::refInfo(const SDictionary& pref) {
	try {
		SeqList reference;
		reference.open(pref["reference"]);
		AnnotDB adb;
		size_t count[3] = { 0, 0, 0 };
		adb.open(pref["annotdb"]);
		auto geneTbl = adb["gene"];
		count[0] = geneTbl.where("type&" + S((int)GENE_TYPE::PROTEIN_CODING)).count(); geneTbl.reset();
		count[1] = geneTbl.where("type&" + S((int)GENE_TYPE::NON_CODING)).count(); geneTbl.reset();
		count[2] = geneTbl.where("type!=0").count();
		if (pref["oformat"] == "_obj_") {
			sobj lgs = SArray();
			sfor(reference) {
				lgs.add({
					D_("name", $_.name),
					D_("length", $_.length()),
					D_("mask", $_.mask.length(true))
					});
			}
			sobj summary = {
				D_("lg", lgs),
				D_("genes", sobj({
					D_("coding", count[0]),
					D_("ncoding", count[1]),
					D_("total", count[2])
					})),
				D_("source", reference.attribute.hasKey("source") ? reference.attribute["source"] : "Unknown"),
				D_("species", reference.attribute.hasKey("species") ? reference.attribute["species"] : "Unknown"),
				D_("version", reference.attribute.hasKey("version") ? reference.attribute["version"] : "Unknown"),
			};
			param.response.attribute["result"] = summary;
		}
		else {
			if (param.output.empty()) param.ostream.setStrOStream(param.response.output);
			else {
				param.ofile.open(param.output, MAKE);
				param.ostream.setFileOStream(param.ofile);
				param.response.output = param.output;
			}
			param.ostream.print(S_(= ) * 60);
			param.ostream.print(sstr::rfill("Linkage group:", ' ', 15), reference.size());
			param.ostream.print(sstr::rfill("Source:", ' ', 15), reference.attribute.hasKey("source") ? reference.attribute["source"] : "Unknown");
			param.ostream.print(sstr::rfill("Species:", ' ', 15), reference.attribute.hasKey("species") ? reference.attribute["species"] : "Unknown");
			param.ostream.print(sstr::rfill("Version:", ' ', 15), reference.attribute.hasKey("version") ? reference.attribute["version"] : "Unknown");
			param.ostream.print(NL);
			//
			param.ostream.print(S_(-) * 40);
			param.ostream.print("|", sstr::bfill("Name", ' ', 12),
				"|", sstr::bfill("Size (bp)", ' ', 12),
				"|", sstr::bfill("Mask (bp)", ' ', 12),
				"|");
			param.ostream.print(S_(-) * 40);
			sfor(reference) {
				param.ostream.print("|", sstr::rfill($_.name, ' ', 12),
					"|", sstr::lfill(String($_.length()), ' ', 12),
					"|", sstr::lfill(String($_.mask.length(true)), ' ', 12),
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
				"|", sstr::fill(String(count[2] - count[0] - count[1]), ' ', 12, DIRECTION::BI),
				"|");
			param.ostream.print(S_(-) * 40);
			param.ostream.print(S_(= ) * 60);
			if (param.response.output) SPrint(param.response.output);
		}
	}
	catch (Exception ex) {
		ex.print();
		param.response = ex;
	}
	return param.response;
}