#include "moirei.h"
Response& Moirei::makeOrthoDB(const SDictionary& pref) {
	try {
		SDataBase db;
		if (!param.output.endWith(".db")) param.output += ".db";
		param.response.output = param.output;
		db.open(param.output);
		//
		db.create("taxon", {
			SColumn("id", { D_("type", "integer"), D_("rule", "key")}),
			SColumn("name", {D_("type", "string")}),
			SColumn("species", {D_("type", "string")})
			});
		db.create("ortholog", {
			SColumn("taxon1", { D_("type", "integer")}),
			SColumn("gid1", {D_("type", "string")}),
			SColumn("gene1", {D_("type", "string")}),
			SColumn("taxon2", { D_("type", "integer")}),
			SColumn("gid2", {D_("type", "string")}),
			SColumn("gene2", {D_("type", "string")}),
			SColumn("flag", {D_("type", "integer")}),
			});
		//
		SArray taxons;
		SFile f(pref["input"]);
		// Read header
		auto& ln = param.ln;
		while (f) {
			f.readLine(ln);
			if (ln.empty()) continue;
			if (ln[0] == '#') {
				if (ln.beginWith("# Taxon")) {
					ln.clip(12);
					ln.replace("NCBITaxon:", "");
					auto vals = ln.split(",");
					sfor(vals) taxons.add({ $_.intValue(), "", "" });
				}
				else if (ln.beginWith("# Species")) {
					ln.clip(10);
					auto vals = ln.split(",");
					sfor2(vals, taxons) {
						if ($_1 == "Homo sapiens") $_2[1] = "human";
						else if ($_1 == "Rattus norvegicus") $_2[1] = "rat";
						else if ($_1 == "Mus musculus") $_2[1] = "mouse";
						else if ($_1 == "Danio rerio") $_2[1] = "zebrafish";
						else if ($_1 == "Xenopus tropicalis") $_2[1] = "tfrog";
						else if ($_1 == "Xenopus laevis") $_2[1] = "afrog";
						else if ($_1 == "Drosophila melanogaster") $_2[1] = "fly";
						else if ($_1 == "Caenorhabditis elegans") $_2[1] = "worm";
						else if ($_1 == "Saccharomyces cerevisiae") $_2[1] = "yeast";
						else $_2[1] = "-";
						$_2[2] = $_1;
					}
				}
				continue;
			}
			else break;
		}
		auto taxonTable = db["taxon"];
		taxonTable.insertAll(taxons);
		//
		SArray orthologs;
		// Read body
		while (f) {
			f.readLine(ln);
			if (ln.empty()) continue;
			ln.replace("NCBITaxon:", "");
			ln.replace("Xenbase:", "");
			ln.replace("ZFIN:", "");
			ln.replace("FB:", "");
			ln.replace("WB:", "");
			auto vals = ln.split(TAB);
			if (vals.size() != 13) continue;
			orthologs.add({
				vals[2].intValue(), vals[0], vals[1], vals[6].intValue(), vals[4], vals[5], (vals[11] == "Yes" ? 0x01 : 0x00)
				});
		}
		auto orthoTable = db["ortholog"];
		orthoTable.insertAll(orthologs);
	}
	catch (Exception ex) {
		ex.print();
		return param.response = ex;
	}
	return param.response;
}
Response& Moirei::makeDiseaseDB(const SDictionary& pref) {
	try {
		SDataBase db;
		if (!param.output.endWith(".db")) param.output += ".db";
		param.response.output = param.output;
		db.open(param.output);
		//
		db.create("disease", {
			SColumn("id", { D_("type", "string")}),
			SColumn("name", {D_("type", "string")}),
			SColumn("taxon", {D_("type", "integer")}),
			SColumn("gid", {D_("type", "string")}),
			SColumn("evidence", {D_("type", "string")}),
			SColumn("source", {D_("type", "string")}),
			SColumn("reference", {D_("type", "string")})
			});
		//
		SArray disease;
		SFile f(pref["input"]);
		// Read header
		auto& ln = param.ln;
		while (f) {
			f.readLine(ln);
			if (ln.empty() || ln[0] == '#') continue;
			else break;
		}
		// Read body
		while (f) {
			f.readLine(ln);
			if (ln.empty()) continue;
			auto vals = ln.split(TAB);
			if (vals.size() != 18 || vals[2] != "gene") continue;
			vals[0].replace("NCBITaxon:", "");
			vals[3].replace("Xenbase:", "");
			vals[3].replace("ZFIN:", "");
			vals[3].replace("FB:", "");
			vals[3].replace("WB:", "");
			disease.add({
				vals[6], vals[7], vals[0].intValue(), vals[3], vals[14], vals[17], vals[15]
				});
		}
		//
		auto disTable = db["disease"];
		disTable.insertAll(disease);
	}
	catch (Exception ex) {
		ex.print();
		return param.response = ex;
	}
	return param.response;
}
