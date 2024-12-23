#include "moirei.h"
Response& Moirei::genomeSeq(const SDictionary& pref) {
	try {
		param.setPref(pref);
		SeqList reference;
		reference.open(pref["reference"]);
		auto nameindex = reference.nameIndex();
		srange range(1, 0x7FFFFFFF);
		if (pref.hasKey("limit-query-length")) range.end = pref["limit-query-length"].intValue();
		int climit = 0x7FFFFFFF;
		if (pref.hasKey("limit-query-count")) climit = pref["limit-query-count"].intValue();
		if (pref["verbose"]) SPrint("Available query size : ", range.begin, "-", range.end);
		Array<RefPos> sites(pref["_args_"].size());
		sfor2(sites, pref["_args_"].array()) {
			$_1 = RefPos::toPos($_2, nameindex);
			auto len = $_1.length(true);
			if (len < range.begin || range.end < len) throw RangeException(outRangeErrorText("Sequence length", len, range.begin, range.end));
		}
		if (sites.empty() || climit < sites.size()) throw RangeException(outRangeErrorText("Query count", sites.size(), 1, climit));

		SeqList seqs;
		seqs.resize(sites.size());
		// Reference
		if (pref["verbose"]) SPrint("Load reference.");
		reference.load(pref["reference"]);
		// Annotation
		sushort atypes = 0;
		AnnotDB db;
		if (pref.hasKey("annotation") && pref["annotation"].size()) {
			if (pref["verbose"]) SPrint("Connect annotation database.");
			if (!pref.hasKey("annotdb")) throw InsufficientArgsException(insufficientArgsErrTxt("annotdb"));
			db.open(pref["annotdb"]);
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
				case 'M':
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
		if (pref["verbose"]) SPrint("Annotation flag : ", (int)atypes);


		// Extraction
		sfor2(sites, seqs) {
			if ($_1.idx == -1) continue;
			if (pref["verbose"]) SPrint("Get sequence for '", $_1.toString(reference), "'.");
			if (atypes) db.annotate(reference[$_1.idx], $_1, atypes);
			$_2 = reference.subsequence($_1);
			if ($_1.dir) $_2.complement();
		}
		// Output
		if (pref["verbose"]) SPrint("Completed.");
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
