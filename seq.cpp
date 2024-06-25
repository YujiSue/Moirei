#include "sapp/scuiapp.h"
#include "sbioinfo/annotation.h"
#include "moirei.h"
using namespace slib;
using namespace slib::smath;
using namespace slib::sio;
using namespace slib::sutil;
using namespace slib::sbio;
using namespace slib::sbio::sio;
using namespace slib::sbio::sutil;

Response moir::refSearch(moir::Param& param, const SDictionary& pref) {
	Response res;
	SeqSearchParam par(DNA_SEQ4);
	//
	if (pref.hasKey("seed-length")) par.setSeed(pref["seed-length"]);
	if (pref.hasKey("min-match")) par.min_match = pref["min-match"];
	if (pref.hasKey("max-gap")) par.max_gap = pref["max-gap"];
	if (pref.hasKey("max-missmatch")) par.max_miss = pref["max-missmatch"];
	if (pref.hasKey("min-score")) par.min_score = pref["min-score"];
	if (pref.hasKey("ext-threshold")) par.ext_threshold = pref["ext-threshold"];
	//if (pref.hasKey("not-complement-search")) par.complement = false;
	//if (pref.hasKey("strict-align")) par.strict = true;
	if (pref.hasKey("thread")) par.setAsync(pref["thread"]);
	//
	bool exact_match = false;
	if (pref.hasKey("exact-match")) exact_match = pref["exact-match"];
	//
	bool show_align = false;
	if (pref.hasKey("show-alignment")) show_align = pref["show-alignment"];
	//
	SeqList reference;
	reference.load(pref["reference"]);
	auto refnum = reference.size();
	//
	SeqSearch searcher(&par);
	DNASeqTrie2 trie(&par);
	//
	if (pref.hasKey("input")) {
		Fasta fa(pref["input"], sseq::DNA);
		fa.makeIndex();
		Sequence seq;
		trie.resize(fa.count());
		sforin(s, 0, fa.count()) {
			fa >> seq;
			trie.queries[2 * s];
		}
	}
	else {
		sforeach(query, pref["_args_"]) trie.addQuery(query.toString());
	}
	trie.complete();
	searcher.search(reference, trie);
	auto qnum = trie.qcount();
	sforin(q, 0, qnum) {
		auto row = searcher.aligns[q];
		sforin(r, 0, refnum) {
			sforeach(align, row[r]) {
				if ((exact_match && 
					align.cigars[0].option == scigar::PMATCH && 
					align.cigars[0].length == trie.queries[q].size()) ||
					(!exact_match && par.min_score <= align.score)) {
					//
					param.ostream << reference[align.ref.idx].name << ":" <<
						align.ref.begin + 1 << "-" << align.ref.end + 1 <<
						"(" << (align.ref.dir ? "-" : "+") << ")" << "  " <<
						align.cigars.toString() << NL;
					//
					if (show_align) {




					}
				}
			}
		}
	}
	if (res.code) throw SAppException(res);
	return res;
}


Response moir::primerInfo(const SDictionary& pref) {
	Response res;
	SeqSearchParam par(DNA_SEQ4);



	if (res.code) throw SAppException(res);
	return res;
}