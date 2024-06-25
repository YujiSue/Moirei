#include "moirei.h"
using namespace slib;
using namespace slib::sbio;


Response moir::getProtein(moir::Param& param, const SDictionary& pref) {
	Response res;
	SeqList reference(pref["reference"], false), aaseqs;
	Sequence seq(DNA_SEQ);
	AnnotDB db(pref["annotdb"]);
	auto codon = DEFAULT_CODON;
	if (pref["codon-table"]) {
		auto codes = sjson::parse(pref["codon-table"]);
		sforin(i, 0, 4) {
			sforin(j, 0, 4) {
				sforin(k, 0, 4) {
					codon[i][j][k] = codes[i][j][k].intValue();
				}
			}
		}
	}
	if (pref["query-type"] == "gene") {
		/*
		sfor(pref["_args_"]) {
			auto &genes = db.searchGenes($_);
			if (genes.size()) {
				sforeach(ginfo, genes) {
					if (ginfo.type & (int)GENE_TYPE::PROTEIN_CODING) {
						sforeach(tinfo, ginfo.transcripts) {
							if (tinfo->type == M_RNA) {
								seq.name = sstr::toUpper(tinfo->name);
								seq.setSeqAs("", DNA_SEQ);
								sforeach(sinfo, tinfo->structures) {
									if (sinfo.type == CDS) seq << reference.raw(RefPos(tinfo->idx, sinfo.begin - 1, sinfo.end - 1));
								}
								if (ginfo.dir) seq.complement();
								aaseqs.add(seq.transcribe().translate(codon));
							}
						}
					}
				}
			}
		}
		*/
	}
	else if (pref["query-type"] == "transcript") {
		sfor(pref["_args_"]) {
			/*
			auto& genes = db.searchGenes($_);
			if (genes.size()) {
				sforeach(gene, genes) {
					sforeach(transcript, gene.transcripts) {
						if (transcript->type == M_RNA) {
							seq.name = sstr::toUpper(transcript->name);
							seq.setSeqAs("", DNA_SEQ);
							sforeach(sinfo, transcript->structures) {
								if (sinfo.type == CDS) seq << reference.raw(RefPos(gene.idx, sinfo.begin, sinfo.end));
							}
							if (gene.dir) seq.complement();
							aaseqs.add(seq.transcribe().translate(codon));
						}
					}
				}
			}
			*/
		}
	}
	String format = pref["oformat"] == "auto" ? (pref["output"] ? sfs::extension(pref["output"]) : "fa") : pref["oformat"].string();
	if (pref["output"]) {
		param.ofile.open(pref["output"], MAKE);
		param.ostream.setFileOStream(param.ofile);
	}
	sbio::sio::writeSeqs(param.ostream, aaseqs, format);
	return res;
}

