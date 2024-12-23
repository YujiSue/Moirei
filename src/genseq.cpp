#include "moirei.h"
Response& Moirei::geneSeq(const SDictionary& pref) {
	try {
		param.setPref(pref);
		//
		SeqList reference, seqs;
		if (pref["verbose"]) SPrint("Load reference.");
		reference.load(pref["reference"]);
		if (pref["verbose"]) SPrint("Connect annotation db.");
		// Annotation
		sushort atypes = 0;
		AnnotDB db;
		if (pref.hasKey("annotation") && pref["annotation"].size()) {
			if (pref["verbose"]) SPrint("Connect annotation database.");
			if (!pref.hasKey("annotdb")) throw InsufficientArgsException(insufficientArgsErrTxt("annotdb"));
			//
			auto str = sstr::toUpper(pref["annotation"]);
			sforeach(c, str) {
				switch (c) {
				case 'G':
					atypes |= (sushort)ANNOT_CATEGORY::GENE;
					break;
				case 'T':
					atypes |= (sushort)ANNOT_CATEGORY::TRANSCRIPT;
					break;
				case 'V':
					atypes |= (sushort)ANNOT_CATEGORY::MUTATION;
					break;
				case 'F':
					atypes |= (sushort)ANNOT_CATEGORY::FEATURE;
					break;
				default:
					break;
				}
			}
		}
		db.open(pref["annotdb"]);
		sfor(pref["_args_"]) {
			if (pref["verbose"]) SPrint("Search gene(s) '", $_, "'.");
			//
			auto genes = db.getGenes($_);
			sforeach(gene, genes) {
				if (atypes) db.annotate(reference[gene.idx], gene, atypes);
				Sequence seq = reference[gene.idx].subsequence(slib::shift(gene, -1));
				seq.name = gene.name;
				if (gene.dir) seq.complement();
				seqs.add(seq);
			}
		}
		if (pref["verbose"]) SPrint("Completed.");
		// Output
		if (param.oformat == "auto") param.oformat = "fa";
		moir::exportSeqs(param, seqs);
		if (param.response.output && !pref["silent"]) SPrint(param.response.output);
	}
	catch (Exception ex) {
		ex.print();
		param.response = ex;
	}
	return param.response;
}
