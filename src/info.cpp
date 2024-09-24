#include "moirei.h"




Response& Moirei::primerInfo(const SDictionary& pref) {
	try {
		String seq = (const char*)pref["_args_"][0], cseq, parts[2];
		seq.trim();
		int len = seq.size();
		int gc = sna::gcCount(seq);
		cseq = sna::complement(seq);
		parts[0] = seq.substring(0, seq.size() / 2);
		parts[1] = seq.substring(seq.size() / 2);
		float bias = abs((float)sna::gcCount(parts[0]) / parts[0].size() - (float)sna::gcCount(parts[1]) / parts[1].size());
		auto mt = sdna::meltTemp(seq);
		
		int dist[2] = {
			smath::levenshtein(parts[0].cstr(), parts[0].size(), cseq.cstr(), cseq.size()) - (cseq.size() - parts[0].size()),
			smath::levenshtein(parts[1].cstr(), parts[1].size(), cseq.cstr(), cseq.size()) - (cseq.size() - parts[1].size())
		};

		SeqList reference;
		SArray offtarget;
		if (pref.hasKey("reference")) {
			reference.load(pref["reference"]);


			

		}

		if (param.oformat == "_obj_") {
			sobj result = {
				D_("query", seq),
				D_("complement", cseq),
				D_("gc", sfrac(gc, len)),
				D_("melt", mt.first),
				D_("temp", sobj({mt.second[0], mt.second[1], mt.second[2]})),
				D_("bias", bias),
				D_("self", sfrac(sstat::getMin(dist[0], dist[1]), len)),
				D_("off", offtarget)
			};

			SPrint(result);

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


			}
			else {

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


Response& Moirei::motifInfo(const SDictionary& pref) {
	try {
		AnnotDB adb;
		
		intarray protids;
		if (pref["query-type"] == "gene") {

		}
		else if (pref["query-type"] == "transcript") {

		}
		else {

		}

	}
	catch (Exception ex) {
		ex.print();
		param.response = ex;
	}
	return param.response;
}
