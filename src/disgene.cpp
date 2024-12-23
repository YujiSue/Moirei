#include "moirei.h"
/**
* DOID
*   |- name: <disease>
*   |- genes : {}
* 
*  Gene ID
*	|- ortholog : []
*         |- geneid
*         |- taxonid
*   |- reference : []
*/

Response& Moirei::diseaseRelated(const SDictionary& pref) {
	try {
		param.setPref(pref);
		//
		SDataBase db;
		db.open(pref["db"]);
		int taxid;
		auto taxtbl = db["taxon"];
		taxtbl
			.where(sdb::textQuery("name", pref["taxon"]), "or", sdb::textQuery("id", pref["taxon"]))
			.select({ "id" });
		taxid = taxtbl.nrow() ? taxtbl[0][0].intValue() : -1;
		if (taxid < 0) throw NotFoundException(nofoundErrorText(pref["taxon"], "Either human, rat, mouse, trog, afrog, zebrafish, fly, worm, or yeast"));
		//
		bool orthosearch = pref.hasKey("search-ortholog") ? (bool)pref["search-ortholog"] : false;
		//
		sobj result = SDictionary();
		String cond;
		sfor(pref["_args_"]) cond << "name like '%" << $_.toString() << "%' AND ";
		if (cond.size()) cond.resize(cond.size() - 5);
		auto distbl = db["disease"];
		distbl
			.where(cond)
			.select({ "id", "name" });
		sforeach(row, distbl) {
			if (!result.hasKey(row[0])) {
				result.set(row[0], sobj({
					D_("name", row[1]),
					D_("genes", SDictionary())
					}));
			}
		}
		auto dids = result.keyset();
		sfor(dids) {
			SDictionary& genes = result[$_]["genes"];
			distbl.clearAll();
			distbl
				.where(sdb::textQuery("id", $_), "and", "taxon=" + S(taxid))
				.select({ "gid", "reference" });
			sforeach(row, distbl) {
				if (!genes.hasKey(row[0])) genes.set(row[0], { D_("ortho", SDictionary()), D_("reference", SArray()) });
				genes[row[0]]["ortho"].set(row[0], taxid);
				genes[row[0]]["reference"].add(row[1]);
			}
			if (orthosearch) {
				distbl.clearAll();
				distbl
					.join("ortholog", sdb::JOIN::INNER, "ortholog.gid1=disease.gid")
					.where(sdb::textQuery("id", $_), "and", "ortholog.taxon2=" + S(taxid))
					.select({ "ortholog.gid2", "disease.gid", "ortholog.taxon1", "disease.reference" });
				sforeach(row, distbl) {
					if (!genes.hasKey(row[0])) genes.set(row[0], { D_("ortho", SDictionary()), D_("reference", SArray())});
					genes[row[0]]["ortho"].set(row[1], row[2]);
					genes[row[0]]["reference"].add(row[3]);
				}
				distbl.clearAll();
				distbl
					.join("ortholog", sdb::JOIN::INNER, "ortholog.gid2=disease.gid")
					.where(sdb::textQuery("id", $_), "and", "ortholog.taxon1=" + S(taxid))
					.select({ "ortholog.gid1", "disease.gid", "ortholog.taxon2", "disease.reference" });
				sforeach(row, distbl) {
					if (!genes.hasKey(row[0])) genes.set(row[0], { D_("ortho", SDictionary()), D_("reference", SArray()) });
					genes[row[0]]["ortho"].set(row[1], row[2]);
					genes[row[0]]["reference"].add(row[3]);
				}
			}
		}
		if (param.oformat == "_obj_") param.response.attribute["result"] = result;
		else {
			if (param.output.empty()) param.ostream.setStrOStream(param.response.output);
			else {
				param.ofile.open(param.output, MAKE);
				param.ostream.setFileOStream(param.ofile);
				param.response.output = param.output;
			}
			if (param.oformat == "tsv") {
				param.ostream.print("DOID", TAB, "Disease", TAB, "Related gene", TAB, "Reported gene(s)", TAB, "Reference");
				sfor(result) {
					SDictionary& genes = $_.value()["genes"];
					sforeach(gene, genes) {
						SArray& references = gene.value()["reference"];
						references.unique();
						SDictionary& orthologs = gene.value()["ortho"];
						String reported = "";
						if (orthologs.hasKey(gene.key())) reported += gene.key() + ",";
						sforeach(orth, orthologs.keyset()) {
							if (orth == gene.key()) continue;
							reported += orth + " (taxon:" + orthologs[orth].toString() + "),";
						}
						if (reported.size()) reported.resize(reported.size() - 1);
						param.ostream << $_.key() << TAB << $_.value()["name"] << TAB << 
							gene.key() << TAB << reported << TAB << references.toString("csv") << NL;
						param.ostream.flush();
					}
				}
			}
			else {
				sfor(result) {
					param.ostream.print(S("=") * 80);
					param.ostream.print(stxt::wrap($_.value()["name"], 80, 20).replace(0, 20, sstr::rfill($_.key(), ' ', 20)));
					param.ostream.print(S("-") * 80);
					param.ostream.print(SP * 1, "Related gene(s):");
					SDictionary& genes = $_.value()["genes"];
					sforeach(gene, genes) {
						param.ostream.print(S("- ") * 40);
						SArray& references = gene.value()["reference"];
						references.unique();
						SDictionary& orthologs = gene.value()["ortho"];
						String reported = "";
						if (orthologs.hasKey(gene.key())) reported += gene.key() + ", ";
						sforeach(orth, orthologs.keyset()) {
							if (orth == gene.key()) continue;
							reported += orth + " (taxon:" + orthologs[orth].toString() + "), ";
						}
						if (reported.size()) reported.resize(reported.size() - 2);
						reported << " [" << references.toString(", ") << "]";
						param.ostream.print(SP * 3, gene.key());
						param.ostream.print(stxt::wrap(reported, 80, 15).replace(3, 10, "Reference:"));
					}
					param.ostream.print(S("=") * 80);
				}
			}
			if (param.response.output) SPrint(param.response.output);
		}
	}
	catch (Exception ex) {
		ex.print();
		param.response = ex;
	}
	return param.response;
}

Response& Moirei::diseaseInfo(const SDictionary& pref) {
	try {
		param.setPref(pref);
		//
		SDataBase db;
		db.open(pref["db"]);
		//
		sobj result = SDictionary();
		stringarray query;
		sfor(pref["_args_"]) query.add($_.toString());
		auto distbl = db["disease"];
		distbl
			.where(sdb::textsQuery("gid", query, MATCH::PARTIAL, false))
			.select({ "id", "name", "gid", "reference"});
		sforeach(row, distbl) {
			if (!result.hasKey(row[0])) {
				result.set(row[0], sobj({
					D_("name", row[1]),
					D_("genes", SArray()),
					D_("reference", SArray())
					}));
			}
			result[row[0]]["genes"].append(row[2]);
			result[row[0]]["genes"].array().unique();
			result[row[0]]["reference"].append(row[3]);
			result[row[0]]["reference"].array().unique();
		}
		if (param.oformat == "_obj_") param.response.attribute["result"] = result;
		else {
			if (param.output.empty()) param.ostream.setStrOStream(param.response.output);
			else {
				param.ofile.open(param.output, MAKE);
				param.ostream.setFileOStream(param.ofile);
				param.response.output = param.output;
			}
			if (param.oformat == "tsv") {
				param.ostream.print("DOID", TAB, "Disease", TAB, "Gene", TAB, "Reference");
				sfor(result) {
					auto& dict = $_.value().dict();
					param.ostream.print($_.key(), TAB, dict["name"], TAB, dict["genes"].toString("csv"), TAB, dict["reference"].toString("csv"));
				}
			}
			else {
				sfor(result) {
					auto& dict = $_.value().dict();
					param.ostream.print(S("=") * 80);
					param.ostream.print(stxt::wrap($_.value()["name"], 80, 20).replace(0, 20, sstr::rfill($_.key(), ' ', 20)));
					param.ostream.print(S("-") * 80);
					param.ostream.print(SP * 1, "Gene(s):");
					param.ostream.print(S("- ") * 40);
					param.ostream.print(SP * 3, dict["genes"].toString("csv"));
					param.ostream.print(stxt::wrap(dict["reference"].toString(", "), 80, 15).replace(3, 10, "Reference:"));
					param.ostream.print(S("=") * 80);
				}
			}
			if (param.response.output) SPrint(param.response.output);
		}
	}
	catch (Exception ex) {
		ex.print();
		param.response = ex;
	}
	return param.response;
}

Response& Moirei::relatedDisease(const SDictionary& pref) {
	try {
		param.setPref(pref);
		//
		AnnotDB adb;




		SDataBase db;
		db.open(pref["db"]);
		int taxid;
		auto taxtbl = db["taxon"];
		taxtbl
			.where(sdb::textQuery("name", pref["taxon"]), "or", sdb::textQuery("id", pref["taxon"]))
			.select({ "id" });
		taxid = taxtbl.nrow() ? taxtbl[0][0].intValue() : -1;
		if (taxid < 0) throw NotFoundException(nofoundErrorText(pref["taxon"], "Either human, rat, mouse, trog, afrog, zebrafish, fly, worm, or yeast"));
		//
		bool orthosearch = pref["search-ortholog"];
		//
		sobj result = SDictionary();
		String cond;
		sfor(pref["_args_"]) cond << "name like '%" << $_.toString() << "%' AND ";
		if (cond.size()) cond.resize(cond.size() - 5);
		auto distbl = db["disease"];
		distbl
			.where(cond)
			.select({ "id", "name" });
		sforeach(row, distbl) {
			if (!result.hasKey(row[0])) {
				result.set(row[0], sobj({
					D_("name", row[1]),
					D_("genes", SDictionary())
					}));
			}
		}
		auto dids = result.keyset();
		sfor(dids) {
			SDictionary& genes = result[$_]["genes"];
			distbl.clearAll();
			distbl
				.where(sdb::textQuery("id", $_), "and", "taxon=" + S(taxid))
				.select({ "gid", "reference" });
			sforeach(row, distbl) {
				if (!genes.hasKey(row[0])) genes.set(row[0], { D_("ortho", SDictionary()), D_("reference", SArray()) });
				genes[row[0]]["ortho"].set(row[0], taxid);
				genes[row[0]]["reference"].add(row[1]);
			}
			if (orthosearch) {
				distbl.clearAll();
				distbl
					.join("ortholog", sdb::JOIN::INNER, "ortholog.gid1=disease.gid")
					.where(sdb::textQuery("id", $_), "and", "ortholog.taxon2=" + S(taxid))
					.select({ "ortholog.gid2", "disease.gid", "ortholog.taxon1", "disease.reference" });
				sforeach(row, distbl) {
					if (!genes.hasKey(row[0])) genes.set(row[0], { D_("ortho", SDictionary()), D_("reference", SArray()) });
					genes[row[0]]["ortho"].set(row[1], row[2]);
					genes[row[0]]["reference"].add(row[3]);
				}
				distbl.clearAll();
				distbl
					.join("ortholog", sdb::JOIN::INNER, "ortholog.gid2=disease.gid")
					.where(sdb::textQuery("id", $_), "and", "ortholog.taxon1=" + S(taxid))
					.select({ "ortholog.gid1", "disease.gid", "ortholog.taxon2", "disease.reference" });
				sforeach(row, distbl) {
					if (!genes.hasKey(row[0])) genes.set(row[0], { D_("ortho", SDictionary()), D_("reference", SArray()) });
					genes[row[0]]["ortho"].set(row[1], row[2]);
					genes[row[0]]["reference"].add(row[3]);
				}
			}
		}
		if (param.oformat == "_obj_") param.response.attribute["result"] = result;
		else {
			if (param.output.empty()) param.ostream.setStrOStream(param.response.output);
			else {
				param.ofile.open(param.output, MAKE);
				param.ostream.setFileOStream(param.ofile);
				param.response.output = param.output;
			}
			if (param.oformat == "tsv") {
				param.ostream.print("DOID", TAB, "Disease", TAB, "Related gene", TAB, "Reported gene(s)", TAB, "Reference");
				sfor(result) {
					SDictionary& genes = $_.value()["genes"];
					sforeach(gene, genes) {
						SArray& references = gene.value()["reference"];
						references.unique();
						SDictionary& orthologs = gene.value()["ortho"];
						String reported = "";
						if (orthologs.hasKey(gene.key())) reported += gene.key() + ",";
						sforeach(orth, orthologs.keyset()) {
							if (orth == gene.key()) continue;
							reported += orth + " (taxon:" + orthologs[orth].toString() + "),";
						}
						if (reported.size()) reported.resize(reported.size() - 1);
						param.ostream << $_.key() << TAB << $_.value()["name"] << TAB <<
							gene.key() << TAB << reported << TAB << references.toString("csv") << NL;
						param.ostream.flush();
					}
				}
			}
			else {
				sfor(result) {
					param.ostream.print(S("=") * 80);
					param.ostream.print(stxt::wrap($_.value()["name"], 80, 20).replace(0, 20, sstr::rfill($_.key(), ' ', 20)));
					param.ostream.print(S("-") * 80);
					param.ostream.print(SP * 1, "Related gene(s):");
					SDictionary& genes = $_.value()["genes"];
					sforeach(gene, genes) {
						param.ostream.print(S("- ") * 40);
						SArray& references = gene.value()["reference"];
						references.unique();
						SDictionary& orthologs = gene.value()["ortho"];
						String reported = "";
						if (orthologs.hasKey(gene.key())) reported += gene.key() + ", ";
						sforeach(orth, orthologs.keyset()) {
							if (orth == gene.key()) continue;
							reported += orth + " (taxon:" + orthologs[orth].toString() + "), ";
						}
						if (reported.size()) reported.resize(reported.size() - 2);
						reported << " [" << references.toString(", ") << "]";
						param.ostream.print(SP * 3, gene.key());
						param.ostream.print(stxt::wrap(reported, 80, 15).replace(3, 10, "Reference:"));
					}
					param.ostream.print(S("=") * 80);
				}
			}
			if (param.response.output) SPrint(param.response.output);
		}
	}
	catch (Exception ex) {
		ex.print();
		param.response = ex;
	}
	return param.response;
}
