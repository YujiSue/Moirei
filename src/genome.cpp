#include "moirei.h"
using namespace slib;
using namespace slib::sio;
using namespace slib::sutil;
using namespace slib::sbio;
using namespace moir;


Response &Moirei::countGCRatio(const SDictionary& pref) {
	try {
		SeqList reference;
		ubytearray seq;
		if (pref["verbose"]) SPrint("Load reference.");
		reference.load(pref["reference"]);
		if (param.output.empty()) {
			param.outdir = sfs::splitPath(pref["reference"]).first;
			param.output = sfs::joinPath(param.outdir, sfs::fileName(pref["reference"], false) + "_gc.bin");
		}
		param.ofile.open(param.output, MAKE);
		//
		param.ofile.writeInt((int)reference.size());
		int bin = pref["bin"], n = bin / 4;
		bin = n * 4;
		seq.resize(bin);
		//
		param.ofile.writeInt(bin);
		sfori(reference) {
			//
			param.ofile.writeInt(i);
			//
			param.ofile.writeInt((int)((reference[i].size() - 1) / n + 1));
			auto bp = reference[i].data();
			size_t current = 0, end = reference[i].size() - n;
			while (current < end) {
				sdna::expand4(bp, 0, bin, seq.data());
				//
				param.ofile.writeInt((int)sna::gcCount(seq));
				bp += n; current += n;
			}
			seq.resize(reference[i].size() - current);
			sdna::expand4(bp, 0, seq.size(), seq.data());
			//
			param.ofile.writeInt((int)sna::gcCount(seq));
		}
	}
	catch (Exception ex) {
		ex.print();
		param.response = ex;
	}
	return param.response;
}