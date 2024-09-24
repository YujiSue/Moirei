#include "moirei.h"

/**
* 
*/
constexpr subyte ANALYSIS_COMPLETED = 0x01;
constexpr subyte SHORT_SEQUENCE_ERROR = 0x02;
constexpr subyte NOT_ALIGNED_ERROR = 0x04;
constexpr subyte TOO_MUCH_N_ERROR = 0x08;
constexpr subyte NO_VARIANT_ERROR = 0x10;
constexpr subyte NO_SAMPLE_ERROR = 0x20;
constexpr subyte MULTI_HIT_ERROR = 0x40;


typedef int getFilterInfo(slib::SDictionary &info);
typedef int customVarFilter(slib::sbio::Variant *var, sobj opts);
Response& Moirei::varFilter(const SDictionary& pref) {
	SeqList ref;
	ref.load(pref["reference"]);
	//
	AnnotDB adb;
	bool annotation = pref.hasKey("annotdb") && sfs::exist(pref["annotdb"]);
	if (annotation) {
		adb.open(pref["annotdb"]);
		adb.loadGenes();
		adb.loadMutations();
		adb.loadVariants();
	}
	//
	bool novel = pref.hasKey("novel-only") && pref["novel-only"];
	//
	param.oformat = pref.hasKey("oformat") ? (pref["oformat"] == "auto" ? "tsv" : (const char*)pref["oformat"]) : "tsv";
	//
	VarParam vpar;
	vpar.annot = (pref.hasKey("cds-only") && pref["cds-only"]) ? CDS : 0xFF;
	vpar.homo_select = pref.hasKey("homo-only") && pref["homo-only"];
	if (pref.hasKey("min-qual")) {
		auto val = pref["min-qual"].floatValue();
		sforin(i, 0, 3) vpar.vcp.min_qual[i] = val;
		vpar.svp.min_qual = val;
	}
	if (pref.hasKey("min-freq")) {
		auto val = pref["min-freq"].floatValue();
		sforin(i, 0, 3) vpar.vcp.min_freq[i] = val;
		vpar.svp.min_freq = val;
	}
	//
	bool plf = false, plio = false;
	sapp::SPlugIn<slib::sbio::Variant*, sobj> plgf;
	sapp::SPlugIn<const char*, slib::sbio::VarList*, sobj> plgio;
	sobj pfopts, piopts;
	if (pref.hasKey("plugin-filter")) {
		auto pargs = pref["plugin-filter"].split("?");
		try {
			plgf.load(pargs[0]); 
			plgf.call("customVarFilter");
			if (1 < pargs.size() && pargs[1].size()) pfopts = pargs[1].parse("&", "=");
			plf = true;
		}
		catch (Exception ex) { ex.print(); }
	}
	if (pref.hasKey("plugin-io")) {
		auto pargs = pref["plugin-io"].split("?");
		try {
			plgio.load(pargs[0]);
			plgio.call("customVarIO");
			if (1 < pargs.size() && pargs[1].size()) piopts = pargs[1].parse("&", "=");
			plio = true;
		}
		catch (Exception ex) { ex.print(); }
	}
	//
	Variant* var;
	VarFilter filter(&ref, &adb, &vpar);
	//
	VarList vlist;
	sfor(pref["_args_"]) {
		auto path = sfs::splitPath($_);
		vlist.load($_);
		sforeach(var, vlist) {
			if (annotation && (var->varid.empty() || var->varid == ".")) {
				auto verified = adb.verifyVariant(*var, vpar);
				if (novel && verified) var->flag |= NOT_USE_FLAG;
			}
			filter.check(var);
			if (plf) plgf.exec(var, pfopts);
		}
		vlist.tidyUp();
		auto output = sfs::joinPath(path.first, path.second.replace(sfs::extension($_), "filtered." + param.oformat));
		vlist.save(output, "cols=VarID,Chr1,Pos1,Len1,Ref,Alt/Ins,Type,Qual,Freq,Genotype,Gene,Site,Mutation,Substitution,Filter");
		SPrint("Save result to '", output, "'");
		if (plio) { plgio.exec($_, &vlist, piopts); }
		vlist.clearAll();
	}
	return param.response;
}