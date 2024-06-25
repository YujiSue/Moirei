#include "moirei.h"
using namespace slib;
using namespace slib::sio;
using namespace slib::sutil;
using namespace slib::sbio;
using namespace moir;
Response moir::makeDiseaseDB(moir::Param& param, const SDictionary& pref) {
	try {
		String source = pref["source"] ? (const char*)pref["source"] : "";
		Response res;
		SDataBase db(pref["output"]);
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
		SFile f(pref["list"]);
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
		return res;
	}
	catch (Exception ex) {
		ex.print();
		return Response(ex);
	}
}
Response moir::diseaseRelatedGene(moir::Param& param, const SDictionary& pref) {
	try {
		Response res;
		SDataBase db(pref["dbpath"]);
		//SDataBase xdb(sfs::joinPath(pref["refdir"], "ortho.db")), disease(sfs::joinPath(pref["refdir"], "disease.db"));
		SArray specs = pref["species"].split("+");
		sfor(specs) {
			if ($_ == "HS") $_ = "Human";
			if ($_ == "MM") $_ = "Mouse"; 
			if ($_ == "DM") $_ = "Fruit fly";
			if ($_ == "CE") $_ = "Nematode";
			if ($_ == "SC") $_ = "Yeast";
		}
		auto taxons = db["taxon"].select({"*"});

		auto dis = db["disease"];
		//dis.where(sdb::textsQuery("name", )):


		auto queries = pref["disease"].split("+");
		dis.where(sdb::textsQuery("name", queries, MATCH::PARTIAL, true))
			.grouping({ "gid" })
			.select({ "id,name,taxon,gid,evidence,reference" });
		if (pref["search-ortho"]) {


		}
		else {



		}

		return res;
	}
	catch (Exception ex) {
		ex.print();
		return Response(ex);
	}
}