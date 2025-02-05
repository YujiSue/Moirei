#include "moirei.h"
Response& Moirei::geneInfo(const SDictionary& pref) {
	try {
		AnnotDB adb;
		adb.open(pref["annotdb"]);
		auto& genes = adb.getGenes(pref["_args_"][0], MATCH::EXACT, {
			D_("synonym", true),
			D_("xref",  true),
			D_("transcript", pref.hasKey("transcript") ? true : false)
			});
		if (pref["oformat"] == "_obj_") {
			sobj objs = SArray();
			sforeach(gene, genes) {
				objs.add({
					D_("gid", gene.geneid)
					});
			}
			param.response.attribute["result"] = objs;
		}
		else {
			if (param.output.empty()) param.ostream.setStrOStream(param.response.output);
			else {
				param.ofile.open(param.output, MAKE);
				param.ostream.setFileOStream(param.ofile);
				param.response.output = param.output;
			}
			if (param.oformat == "tsv") {
				sforeach(gene, genes) {
					param.ostream.print("Gene ID", TAB, gene.geneid);
					param.ostream.print("Symbol", TAB, gene.name);
					param.ostream.print("Synonyms", TAB, toString(gene.synonym));
					param.ostream.print("Type", TAB, sbio::sutil::geneType(gene.type));
					param.ostream.print("Position",TAB, adb.chromosomes[gene.idx].name, ":", sbio::sutil::gbkPos(gene));
					param.ostream.print("Description", TAB, gene.description);
					param.ostream.print("Attribute");
					sforeach(attr, gene.attribute) {
						param.ostream.print(TAB, attr.key(), TAB, attr.value());
					}
					param.ostream.print("Transcript");
					param.ostream.print(TAB, "#", TAB, "Name", TAB, "Type", TAB, "Size", TAB, "Exon", TAB, "CDS");
					sforeach(rna, gene.transcripts) {
						param.ostream.print(TAB,
							INDEX_(rna, gene.transcripts) + 1, TAB,
							rna->name, TAB,
							sbio::sutil::transcriptType(rna->type), TAB,
							S(rna->length(true)), TAB,
							sbio::sutil::gbkPos(rna->exons(), -1, gene.dir), TAB,
							(rna->type == (int)TRANSCRIPT_TYPE::M_RNA ? sbio::sutil::gbkPos(rna->coding(),-1,  gene.dir) : ""));
					}
				}
			}
			else {
				sforeach(gene, genes) {
					param.ostream.print(S("=") * 90);
					param.ostream.print(sstr::rfill("Gene ID:", ' ', 15), gene.geneid);
					param.ostream.print(sstr::rfill("Symbol:", ' ', 15), gene.name);
					param.ostream.print(sstr::rfill("Synonyms:", ' ', 15), toString(gene.synonym));
					param.ostream.print(sstr::rfill("Type:", ' ', 15), sbio::sutil::geneType(gene.type));
					param.ostream.print(sstr::rfill("Position:", ' ', 15), adb.chromosomes[gene.idx].name, ":", gene.begin, "-", gene.end, " (", (gene.dir?"-":"+"), ")");
					param.ostream.print(stxt::wrap(gene.description, 90, 15).replace(0, 12, "Description:"));
					param.ostream.print("Attribute:");
					param.ostream.print(SP * 15, S("-") * 63);
					param.ostream.print(SP * 15, "|", sstr::bfill("KEY", ' ', 15), "|", sstr::bfill("VALUE", ' ', 45), "|");
					param.ostream.print(SP * 15, "|", S("-") * 15, "|", S("-") * 45, "|");
					sforeach(attr, gene.attribute) {
						param.ostream.print(SP * 15, "| ", sstr::rfill(attr.key(), ' ', 14), "| ", sstr::rfill(attr.value(), ' ', 44), "|");
					}
					param.ostream.print(SP * 15, S("-") * 63);
					param.ostream.print("Transcript:");
					param.ostream.print(SP * 15, S("-") * 69);
					param.ostream.print(SP * 15, "|", sstr::bfill("Name", ' ', 26), "|", sstr::bfill("Type", ' ', 21), "|", sstr::bfill("Size (bp)", ' ', 18), "|");
					param.ostream.print(SP * 15, "|", S("-") * 26, "|", S("-") * 21, "|", S("-") * 18, "|");
					sforeach(rna, gene.transcripts) {
						param.ostream.print(SP * 15, "| ", 
							sstr::rfill(rna->name, ' ', 25), "| ", 
							sstr::rfill(sbio::sutil::transcriptType(rna->type), ' ', 20), "| ", 
							sstr::rfill(S(rna->length(true)), ' ', 17), "|");
					}
					param.ostream.print(SP * 15, S("-") * 69);
					param.ostream.print(S("=") * 90);
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
Response& Moirei::orthoGene(const SDictionary& pref) {
	try {
		param.setPref(pref);
		//
		SDataBase db;
		db.open(pref["db"]);
		int taxid;
		sindex taxlist;
		auto taxtbl = db["taxon"];
		taxtbl.select({ "id", "name" });
		sfor(taxtbl) taxlist[$_["name"]] = $_["id"].intValue();
		if (sstr::isNumeric(pref["taxon"])) taxid = pref["taxon"].intValue();
		else taxid = taxlist[pref["taxon"]];
		if (taxlist.rlookup(taxid).empty()) throw NotFoundException(nofoundErrorText(pref["taxon"], "Either human, rat, mouse, trog, afrog, zebrafish, fly, worm, or yeast"));
		SArray genes;
		String query = (const char*)pref["_args_"][0];
		//
		auto orthtbl = db["ortholog"];
		orthtbl.where("taxon2=" + S(taxid), "AND", "(", sdb::textQuery("gid1", query), "OR", sdb::textQuery("gene1", query, MATCH::PARTIAL), ")")
			.select({ "gid2,gene2,gid1,gene1,taxon1" });
		genes.append(orthtbl.rows());
		orthtbl.clearAll();
		orthtbl.where("taxon1=" + S(taxid), "AND", "(", sdb::textQuery("gid2", query), "OR", sdb::textQuery("gene2", query, MATCH::PARTIAL), ")")
			.select({ "gid1,gene1,gid2,gene2,taxon2" });
		genes.append(orthtbl.rows());
		sobj result = SDictionary();
		sfor(genes) {
			if (!result.hasKey($_[0]))
				result[$_[0]] = { D_("name", $_[1]), D_("origin", SArray()) };
			result[$_[0]]["origin"].add({
				D_("gid", $_[2]),
				D_("gene", $_[3]),
				D_("taxon", taxlist.rlookup($_[4].intValue()))
				});
		}
		//
		if (pref["oformat"] == "_obj_") {
			param.response.attribute["result"] = result;
		}
		else {
			if (param.output.empty()) param.ostream.setStrOStream(param.response.output);
			else {
				param.ofile.open(param.output, MAKE);
				param.ostream.setFileOStream(param.ofile);
				param.response.output = param.output;
			}
			if (param.oformat == "tsv") {
				param.ostream.print("Orthologue", TAB * 2, "Origin", TAB * 2);
				param.ostream.print("Gene ID", TAB, "Symbol", TAB, "Gene ID", TAB, "Symbol");
				sfor(genes) {
					param.ostream.print($_[0], TAB, $_[1], TAB, $_[2], TAB, $_[3]);
				}
			}
			else {
				param.ostream.print(S("=") * 85);
				param.ostream.print("|", sstr::bfill("Orthologue", ' ', 41), "|", sstr::bfill("Origin", ' ', 41), "|");
				param.ostream.print(S(" -") * 42);
				param.ostream.print("|", sstr::bfill("Gene ID", ' ', 20), "|", sstr::bfill("Symbol", ' ', 20),
					"|", sstr::bfill("Gene ID", ' ', 20), "|", sstr::bfill("Symbol", ' ', 20), "|");
				param.ostream.print(S("-") * 85);
				sfor(genes) {
					param.ostream <<
						"|" << ($_[0].size() < 20 ? sstr::bfill($_[0], ' ', 20) : $_[0].toString().substring(0, 17) + "...") <<
						"|" << ($_[1].size() < 20 ? sstr::bfill($_[1], ' ', 20) : $_[1].toString().substring(0, 17) + "...") <<
						"|" << ($_[2].size() < 20 ? sstr::bfill($_[2], ' ', 20) : $_[2].toString().substring(0, 17) + "...") <<
						"|" << ($_[3].size() < 20 ? sstr::bfill($_[3], ' ', 20) : $_[3].toString().substring(0, 17) + "...") <<
						"|" << NL;
				}
				param.ostream.print(S("=") * 85);
			}
			if (param.response.output && !pref["silent"]) SPrint(param.response.output);
		}
	}
	catch (Exception ex) {
		ex.print();
		param.response = ex;
	}
	return param.response;
}
