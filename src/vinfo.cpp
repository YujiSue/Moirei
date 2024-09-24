#include "moirei.h"

Response& Moirei::variantInfo(const SDictionary& pref) {
	try {
		AnnotDB adb;
		adb.open(pref["annotdb"]);
		intarray vars;
		if (pref.hasKey("query-type")) {

		}
		else {
			//auto &vars = adb.

		}
		if (pref["oformat"] == "_obj_") {
			sobj objs = SArray();
			sforeach(vid, vars) {
				//auto vinfo = adb.va

				objs.add({});
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
				
			}
			else if (param.oformat == "vcf") {
				
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