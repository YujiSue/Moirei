#include "moirei.h"
using namespace slib;
using namespace slib::sio;
using namespace slib::sutil;
using namespace slib::sbio;
using namespace moir;
Response moir::geneInfo(moir::Param& param, const SDictionary& pref) {
	try {
		Response res;
		AnnotDB adb(pref["annotdb"]);

		

		return res;
	}
	catch (Exception ex) {
		ex.print();
		return Response(ex);
	}
}
Response moir::makeOrthoDB(moir::Param& param, const SDictionary& pref) {
	try {
		String source = pref["source"] ? (const char*)pref["source"] : "";
		Response res;
		//
		SDataBase db(pref["output"]);
		db.create("taxon", {
			SColumn("id", { D_("type", "integer"), D_("rule", "key")}),
			SColumn("name", {D_("type", "string")}),
			SColumn("species", {D_("type", "string")})
			});
		db.create("orthology", {
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
		SFile f(pref["list"]);
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
					sfor2(vals, taxons) $_2[2] = $_1;
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
		auto orthoTable = db["orthology"];
		orthoTable.insertAll(orthologs);
		return res;
	}
	catch (Exception ex) {
		ex.print();
		return Response(ex);
	}
}
Response moir::getOrthoGenes(moir::Param& param, const SDictionary& pref) {
	try {
		Response res;



		//
		return res;
	}
	catch (Exception ex) {
		ex.print();
		return Response(ex);
	}
}
