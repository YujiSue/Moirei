#include "moirei.h"
Response& Moirei::refSearch(const SDictionary& pref) {
	try {
		SeqSearchParam par(DNA_SEQ4);
		//
		par.setSeed(pref.hasKey("seed-length") ? pref["seed-length"].intValue() : 12);
		par.min_match = pref.hasKey("min-match") ? pref["min-match"].intValue() : 15;
		par.max_gap = pref.hasKey("max-gap") ? pref["max-gap"].intValue() : 2;
		par.max_miss = pref.hasKey("max-missmatch") ? pref["max-missmatch"].intValue() : 2;
		par.min_score = pref.hasKey("min-score") ? pref["min-score"].intValue() : par.min_match * par.apar.pm_score;
		par.ext_threshold = pref.hasKey("threshold") ? pref["threshold"].floatValue() : 1.0;
		par.setAsync(pref.hasKey("thread") ? pref["thread"].intValue() : 1);
		//
		int max_display = pref.hasKey("max-result") ? pref["max-result"].intValue() : 20;
		//
		bool exact_match = false;
		if (pref.hasKey("exact-match")) exact_match = pref["exact-match"];
		//
		SeqList reference;
		reference.load(pref["reference"]);
		auto refnum = reference.size();
		//
		SeqSearch searcher(&par);
		DNASeqTrie2 trie(&par);
		//
		if (pref.hasKey("input")) {
			auto ext = sfs::extension(pref["input"]);
			Sequence seq;
			if (ext == "fa") {
				Fasta fa(pref["input"], sseq::DNA);
				fa.makeIndex();
				sforin(s, 0, fa.count()) {
					fa >> seq;
					trie.addQuery(seq.raw());
				}
			}
			else {
				SFile f(pref["input"]);
				f >> param.ln;
				seq.setSeqAs(param.ln, DNA_SEQ);
				trie.addQuery(seq.raw());
			}
		}
		else {
			sforeach(query, pref["_args_"]) {
				trie.addQuery(query.toString());
			}
		}
		//
		trie.complete();
		searcher.search(reference, trie);
		//
		Array<Array<AlignPair*>> results(trie.qcount() / 2);
		//
		String alignedseq[3];
		//
		sfori(results) {
			//param.ostream.print("> Query #", (i + 1));
			// Sort by score
			sforin(j, 0, 2) {
				auto row = searcher.aligns[2 * i + j];
				sforin(r, 0, refnum) {
					sforeach(align, row[r]) {
						if ((exact_match &&
							align.cigars[0].option == scigar::PMATCH &&
							align.cigars[0].length == trie.queries[i].size()) ||
							(!exact_match && par.min_score <= align.score)) results[i].add(&align);
					}
				}
			}
			results[i].sort([](const AlignPair* a1, const AlignPair* a2) {
				if (a1->score == a2->score) return a1->ref < a2->ref;
				else return a2->score < a1->score;
				});
		}
		// Output
		if (param.oformat == "auto") param.oformat = "txt";
		if (param.oformat == "json") {
			sobj objs = SArray();
			sforin(q, 0, results.size()) {
				param.ln.resize(trie.queries[2 * q].size());
				sdna::decode(trie.queries[2 * q].data(), 0, trie.queries[2 * q].size(), (subyte*)&param.ln[0]);
				sobj matched = SArray();
				sforeach(align, results[q]) {
					alignedseq[0] = align->alref(reference.raw(align->ref));
					alignedseq[1] = align->match();
					alignedseq[2] = align->alque(param.ln.substring(align->query.begin, align->query.length(true)));
					matched.add({
						D_("ref", sobj({
							D_("idx", align->ref.idx),
							D_("name", reference[align->ref.idx].name),
							D_("pos", sobj({ align->ref.begin, align->ref.end})),
							D_("dir", align->ref.dir)
							})),
						D_("pos", sobj({ align->query.begin, align->query.end})),
						D_("pattern", align->cigars.toString()),
						D_("seq", sobj({alignedseq[0], alignedseq[1], alignedseq[2]}))
						});
				}
				objs.add(matched);
			}
			if (param.output.empty()) param.response.attribute["result"] = objs;
			else {
				SPrint("Save to '", param.output, "'.");
				sjson::save(objs, param.output);
			}
		}
		else {
			if (param.output.empty()) param.ostream.setStrOStream(param.response.output);
			else {
				SPrint("Save to '", param.output, "'.");
				param.ofile.open(param.output, MAKE);
				param.ostream.setFileOStream(param.ofile);
			}
			//
			auto qnum = trie.qcount() / 2;
			sforin(q, 0, qnum) {
				param.ostream.print("> Query #", (q + 1));
				// Sort by score
				Array<AlignPair*>& result = results[q];
				// 
				auto count = sstat::getMin(max_display, (int)result.size());
				sforin(r, 0, count) {
					auto& align = *result[r];
					param.ostream.print(SP * 5, sstr::rfill("Reference", ' ', 12),
						reference[align.ref.idx].name, ":", align.ref.begin + 1, "-", align.ref.end + 1, " (", (align.ref.dir ? "-" : "+"), ")");
					param.ostream.print(SP * 5, sstr::rfill("Query", ' ', 12),
						align.query.begin + 1, "-", align.query.end + 1, SP * 2, align.cigars.toString());
					//
					alignedseq[0] = align.alref(reference.raw(align.ref));
					alignedseq[1] = align.match();
					param.ln.resize(align.query.length(true));
					sdna::decode(trie.queries[2 * q].data(), align.query.begin, align.query.length(true), (subyte*)&param.ln[0]);
					alignedseq[2] = align.alque(param.ln);
					//
					param.ostream.print("");
					param.ostream.print(sstr::lfill(S((align.ref.dir ? align.ref.end : align.ref.begin) + 1), ' ', 16), SP, alignedseq[0], SP, S((align.ref.dir ? align.ref.begin : align.ref.end) + 1));
					param.ostream.print(SP * 17, alignedseq[1]);
					param.ostream.print(sstr::lfill(S(align.query.begin + 1), ' ', 16), SP, alignedseq[2], SP, S(align.query.end + 1));
					param.ostream.print(NL);
				}
				if (q < qnum - 1) param.ostream.print(S("-") * 60);
			}
			if (param.response.output.size()) SPrint(param.response.output);
		}
	}
	catch (Exception ex) {
		ex.print();
		param.response = ex;
	}
	return param.response;
}


struct VAPair {
	Array<salign> aligns;
	Array<Pair<int, Variant>> svariants;
	Array<Pair<srange, Variant>> lvariants;
	VAPair() {}
	~VAPair() {}
	float score() const {
		float s = 0.f;
		sfor(aligns) s += $_.score;
		return s;
	}
};
inline void append2map(salign* align, Array<Pair<srange, salign*>>& map, int seqlen) {
	srange range = align->query;
	if (align->ref.dir) sbio::sutil::reverse(range, seqlen);
	if (map.empty() || map[-1].first.end < range.begin) map.add(Pair<srange, salign*>(range, align));
	else if (range.end < map[0].first.begin) map.insert(0, Pair<srange, salign*>(range, align));
	else {
		sfor(map) {
			if ($_.first.overlap(range)) {
				if ($_.first.include(range.begin)) map.insert((int)$INDEX(map) + 1, Pair<srange, salign*>(range, align));
				else map.insert((int)$INDEX(map), Pair<srange, salign*>(range, align));
				break;
			}
			else if (it < map.end() - 1 && it->first.end < range.begin && range.end < (it + 1)->first.begin) {
				map.insert((int)$INDEX(map) + 1, Pair<srange, salign*>(range, align));
				break;
			}
		}
	}
}
inline int appendSize(salign* align, Array<Pair<srange, salign*>>& map, int seqlen) {
	srange range = align->query;
	if (align->ref.dir) sbio::sutil::reverse(range, seqlen);
	if (map.empty() ||
		range.end < map[0].first.begin ||
		map[-1].first.end < range.begin) return align->query.length(true);
	else {
		sfor(map) {
			if (range.end < $_.first.begin) break;
			if ($_.first.include(range)) return 0;
			else if ($_.first.overlap(range)) {
				if ($_.first.include(range.begin)) range.begin = $_.first.end + 1;
				else range.end = $_.first.begin - 1;
			}
		}
	}
	return range.length(true);
}
inline void makeVar(Cigar& cigar, salign* align, Sequence& seq, SeqList& reference, int& rpos, int& qpos, Variant& var) {
	var.flag = SMALL_VARIANT;
	if (cigar.option == scigar::MMATCH) {
		var.type = cigar.length == 1 ? SNV : MNV;
		var.pos[0].idx = align->ref.idx;
		if (align->ref.dir) var.pos[0].begin = rpos - cigar.length + 2;
		else var.pos[0].begin = rpos + 1;
		var.pos[0].end = var.pos[0].begin + cigar.length - 1;
		var.alt = seq.raw(qpos, cigar.length, align->ref.dir);
		var.attribute["_ref_"] = reference[align->ref.idx].raw(var.pos[0].begin - 1, var.pos[0].length(true));
		var.qual = 0.f;
		if (seq.attribute.hasKey("PCON")) {
			auto qual = seq.attribute["PCON"][0]["data"].data();
			sforin(i, 0, cigar.length) var.qual += (float)qual[qpos + i];
			var.qual /= cigar.length;
		}
		else var.qual = 60.f;
		if (align->ref.dir) rpos -= cigar.length;
		else rpos += cigar.length;
		qpos += cigar.length;
	}
	else if (cigar.option == scigar::MATCH) {
		ubytearray tmp(cigar.length);
		sdna::expand4(reference[align->ref.idx].data(), (align->ref.dir ? rpos - cigar.length + 1 : rpos), cigar.length, tmp.data());
		if (align->ref.dir) sna::toComplement16(tmp);
		sforin(l, 0, cigar.length) {
			if (seq[qpos + l] & tmp[l]) seq[qpos + l] -= tmp[l];
		}
		var.genotype = HETERO_VAR;
		var.type = cigar.length == 1 ? SNV : MNV;
		var.pos[0].idx = align->ref.idx;
		if (align->ref.dir) var.pos[0].begin = rpos - cigar.length + 2;
		else var.pos[0].begin = rpos + 1;
		var.pos[0].end = var.pos[0].begin + cigar.length - 1;
		var.alt = seq.raw(qpos, cigar.length, align->ref.dir);
		var.attribute["_ref_"] = reference[align->ref.idx].raw(var.pos[0].begin - 1, var.pos[0].length(true));
		var.qual = 0.f;
		if (seq.attribute.hasKey("PCON")) {
			auto qual = seq.attribute["PCON"][0]["data"].data();
			sforin(i, 0, cigar.length) var.qual += (float)qual[qpos + i];
			var.qual /= cigar.length;
		}
		else var.qual = 60.f;
		if (align->ref.dir) rpos -= cigar.length;
		else rpos += cigar.length;
		qpos += cigar.length;
	}
	else if (cigar.option == scigar::DELETION) {
		var.type = DELETION;
		var.pos[0].idx = align->ref.idx;
		if (align->ref.dir) var.pos[0].begin = rpos - cigar.length + 2;
		else var.pos[0].begin = rpos + 1;
		var.pos[0].end = var.pos[0].begin + cigar.length - 1;
		var.attribute["ref"] = reference[align->ref.idx].raw(var.pos[0].begin - 1, var.pos[0].length(true));
		var.qual = 0.f;
		if (seq.attribute.hasKey("PCON")) {
			auto qual = seq.attribute["PCON"][0]["data"].data();
			var.qual = ((float)qual[qpos - 1] + (float)qual[qpos]) / 2.f;
		}
		else var.qual = 60.f;
		if (align->ref.dir) rpos -= cigar.length;
		else rpos += cigar.length;
	}
	else if (cigar.option == scigar::INSERTION) {
		var.type = INSERTION;
		var.pos[0].idx = align->ref.idx;
		if (align->ref.dir) var.pos[0].begin = rpos - cigar.length + 2;
		else var.pos[0].begin = rpos + 1;
		var.pos[0].end = var.pos[0].begin;
		var.alt = seq.raw(qpos, cigar.length, align->ref.dir);
		var.attribute["ref"] = reference[align->ref.idx].raw(var.pos[0].begin - 2, 1);
		if (var.alt == var.attribute["ref"].string()) var.type = DUPLICATION;
		var.qual = 0.f;
		if (seq.attribute.hasKey("PCON")) {
			auto qual = seq.attribute["PCON"][0]["data"].data();
			var.qual = ((float)qual[qpos - 1] + (float)qual[qpos]) / 2.f;
		}
		else var.qual = 60.f;
		qpos += cigar.length;
	}
}
inline void makeSVar(salign& a1, salign& a2, Variant& var) {
	var.flag = SR_VARIANT;
	if (a1.ref.idx == a2.ref.idx) {
		if (a1.ref.dir == a2.ref.dir) {
			if (a1.ref.dir) {
				if (var.alt.length()) sna::toComplement(var.alt);
				if (a1.ref.begin < a2.ref.end + 1) {
					var.type |= sbio::DUPLICATION;
					var.pos[0] = a1.ref;
					var.pos[0].begin++;
					var.pos[0].end = a2.ref.end + 1;
				}
				else if (a2.ref.end + 1 < a1.ref.begin) {
					var.type |= sbio::DELETION;
					var.pos[0] = a2.ref;
					var.pos[0].begin = a2.ref.end + 2;
					var.pos[0].end = a1.ref.begin;
				}
			}
			else {
				if (a2.ref.begin < a1.ref.begin) {
					var.type |= sbio::DUPLICATION;
					var.pos[0] = a2.ref;
					var.pos[0].begin++;
					var.pos[0].end = a1.ref.end + 1;
				}
				else if (a1.ref.end + 1 < a2.ref.begin) {
					var.type |= sbio::DELETION;
					var.pos[0] = a1.ref;
					var.pos[0].begin = a1.ref.end + 2;
					var.pos[0].end = a2.ref.begin;
				}
			}
		}
		else {
			var.type |= sbio::INVERSION;
			if (a1.ref.dir) {
				var.pos[0] = a1.ref;
				var.pos[0].begin++;
				var.pos[1] = a2.ref;
				var.pos[1].begin++;
			}
			else {
				var.pos[0] = a1.ref;
				var.pos[0].begin = a1.ref.end + 1;
				var.pos[1] = a2.ref;
				var.pos[1].begin = a2.ref.end + 1;
			}
		}
	}
	else {
		var.type |= sbio::TRANSLOCATION;
		if (a1.ref.dir == a2.ref.dir) {
			if (a1.ref.dir) {
				if (var.alt.length()) sna::toComplement(var.alt);
				var.pos[0] = a1.ref;
				var.pos[0].begin++;
				var.pos[1] = a2.ref;
				var.pos[1].begin = a2.ref.end + 1;
			}
			else {
				var.pos[0] = a1.ref;
				var.pos[0].begin = a1.ref.end + 1;
				var.pos[1] = a2.ref;
				var.pos[1].begin++;
			}
		}
		else {
			var.type |= INVERSION;
			if (a1.ref.dir) {
				var.pos[0] = a1.ref;
				var.pos[0].begin++;
				var.pos[1] = a2.ref;
				var.pos[1].begin++;
			}
			else {
				var.pos[0] = a1.ref;
				var.pos[0].begin = a1.ref.end + 1;
				var.pos[1] = a2.ref;
				var.pos[1].begin = a2.ref.end + 1;
			}
		}
	}
}
inline void makeCSVar() {


}
inline void map2res(Array<Pair<srange, salign*>>& map, VAPair& vap) {
	if (map.size() == 1) {
		vap.aligns.add(*map[0].second);
		return;
	}
	sforin(pair, map.begin(), map.end() - 1) {
		auto nxt = pair + 1;
		if (pair->first.overlap(nxt->first)) {
			int len = pair->first.end - nxt->first.begin + 1;
			if (pair->second->ref.dir) {
				sfor(pair->second->cigars) {
					if ($_.option == scigar::MMATCH) {
						len -= $_.length;
						pair->second->ref.begin += $_.length;
						pair->second->query.begin += $_.length;
					}
					else if ($_.option == scigar::DELETION) {
						pair->second->ref.begin += $_.length;
					}
					else if ($_.option == scigar::INSERTION) {
						len -= $_.length;
						pair->second->query.begin += $_.length;
					}
					else {
						if (len <= 0) break;
						else if (len < $_.length) {
							$_.length -= len;
							pair->second->ref.begin += len;
							pair->second->query.begin += len;
						}
						else {
							len -= $_.length;
							pair->second->ref.begin += $_.length;
							pair->second->query.begin += $_.length;
						}
					}
				}
			}
			else {
				srfor(pair->second->cigars) {
					if ($_.option == scigar::MMATCH) {
						len -= $_.length;
						pair->second->ref.end -= $_.length;
						pair->second->query.end -= $_.length;
					}
					else if ($_.option == scigar::DELETION) {
						pair->second->ref.end -= $_.length;
					}
					else if ($_.option == scigar::INSERTION) {
						len -= $_.length;
						pair->second->query.end -= $_.length;
					}
					else {
						if (len <= 0) break;
						else if (len < $_.length) {
							$_.length -= len;
							pair->second->ref.end -= len;
							pair->second->query.end -= len;
						}
						else {
							len -= $_.length;
							pair->second->ref.end -= $_.length;
							pair->second->query.end -= $_.length;
						}
					}
				}
			}
		}
		vap.aligns.add(*pair->second);
	}
	vap.aligns.add(*map[-1].second);
}


inline void fillmap(Array<Pair<srange, salign*>>& map,
	Range<ArrayIterator<AlignPair*>> range,
	Array<Array<Pair<srange, salign*>>>& maps,
	int length,
	SeqSearchParam* par) {
	Array<ArrayIterator<AlignPair*>> nexts;
	size_t max_len = 0, tmp = 0;
	sforin(it, range.begin, range.end) {
		if ($_->ref.length(true) < max_len) break;
		//
		if (par->min_match <= (tmp = appendSize($_, map, length))) {
			if (max_len < tmp) {
				nexts.clear();
				nexts.add($);
				max_len = tmp;
			}
			else if (max_len == tmp) nexts.add($);
		}
	}
	if (nexts.empty()) maps.add(map);
	else if (nexts.size() == 1) {
		append2map(*nexts[0], map, length);
		fillmap(map, Range<ArrayIterator<AlignPair*>>(nexts[0] + 1, range.end), maps, length, par);
	}
	else {
		sfor(nexts) {
			Array<Pair<srange, salign*>> newmap = map;
			append2map(*$_, newmap, length);
			fillmap(map, Range<ArrayIterator<AlignPair*>>($_ + 1, range.end), maps, length, par);
		}
	}
}

inline void findSmallVar(int i, salign& al, Sequence& seq, SeqList& reference, Array<Pair<int, Variant>>& variants) {
	int rpos = al.ref.dir ? al.ref.end : al.ref.begin, qpos = al.query.begin;
	sfor(al.cigars) {
		if ($_.option == scigar::MATCH ||
			$_.option == scigar::MMATCH ||
			$_.option == scigar::DELETION ||
			$_.option == scigar::INSERTION) {
			Variant var;
			makeVar($_, &al, seq, reference, rpos, qpos, var);
			variants.add(Pair<int, Variant>(i, var));
		}
		else {
			if (al.ref.dir) rpos -= $_.length;
			else rpos += $_.length;
			qpos += $_.length;
		}
	}
}
inline void findLargeVariant(int i1, int i2, salign& al1, salign& al2, ubytearray& que, Array<Pair<srange, Variant>>& variants) {
	Variant var;
	srange rng1 = al1.query, rng2 = al2.query;
	//if (al1.ref.dir) sbio::sutil::reverse(rng1, que.size());
	//if (al2.ref.dir) sbio::sutil::reverse(rng2, que.size());
	if (rng1.end + 1 < rng2.begin) {
		var.alt.resize(rng2.begin - rng1.end - 1);
		sdna::decode(que.data(), rng1.end + 1, rng2.begin - rng1.end - 1, (subyte*)&var.alt[0]);
	}
	makeSVar(al1, al2, var);
	variants.add(Pair<srange, Variant>(srange(i1, i2), var));
}
inline void findComplexVariant(Array<Pair<srange, Variant>>& variants) {
	sfor(variants) {
		if (!$_.second.type) continue;
		sforin(vit, $ + 1, variants.end()) {
			if (!vit->second.type) continue;
			Variant var;
			makeCSVar();
			//if (var.type != 0) 

		}
		/*
		for (auto vit_ = it + 1; vit_ < variants.end(); ++vit_) {

			if (((it->var.type & sbio::DELETION && vit_->var.type & sbio::DUPLICATION) ||
				(it->var.type & sbio::DUPLICATION && vit_->var.type & sbio::DELETION)) &&
				it->var.pos[0].idx == vit_->var.pos[0].idx &&
				it->var.pos[0].dir == vit_->var.pos[0].dir) {
				if (it->var.pos[0].dir) {
					if (vit_->var.pos[0].begin <= it->var.pos[0].begin &&
						vit_->var.pos[0].begin + vit_->var.pos[0].end <= it->var.pos[0].begin + it->var.pos[0].end &&
						it->var.pos[0].begin < vit_->var.pos[0].begin + vit_->var.pos[0].end) {
						svar_data cvar;
						int dlen, ilen;
						if (it->var.type & sbio::DELETION) {
							dlen = (it->var.pos[0].begin + it->var.pos[0].end) - (vit_->var.pos[0].begin + vit_->var.pos[0].end);
							ilen = it->var.pos[0].begin - vit_->var.pos[0].begin;
							cvar.type = (0 < dlen ? sbio::DELETION : 0) | sbio::INSERTION;
							cvar.pos[0] = vit_->var.pos[0];
							cvar.pos[0].end = cvar.pos[0].begin + dlen - 1;
							cvar.pos[1] = vit_->var.pos[0];
							cvar.pos[1].end = cvar.pos[1].begin + ilen - 1;
						}
						else {
							dlen = it->var.pos[0].begin - vit_->var.pos[0].begin;
							ilen = (it->var.pos[0].begin + it->var.pos[0].end) - (vit_->var.pos[0].begin + vit_->var.pos[0].end);
							cvar.type = (0 < dlen ? sbio::DELETION : 0) | sbio::INSERTION;
							cvar.pos[0] = vit_->var.pos[0];
							cvar.pos[0].end = cvar.pos[0].begin + dlen - 1;
							cvar.pos[1] = vit_->var.pos[0];
							cvar.pos[1].begin = vit_->var.pos[0].end;
							cvar.pos[1].end = cvar.pos[1].begin + ilen - 1;
						}
						cvar.alt = vit_->var.alt + "/" + it->var.alt;
						it->var = cvar;
						it->idx2 = vit_->idx2;
						vit_->var.type = 0;
					}
				}
				else {
					if (it->var.pos[0].begin <= vit_->var.pos[0].begin &&
						it->var.pos[0].begin + it->var.pos[0].end <= vit_->var.pos[0].begin + vit_->var.pos[0].end &&
						vit_->var.pos[0].begin < it->var.pos[0].begin + it->var.pos[0].end) {
						svar_data cvar;
						int dlen, ilen;
						if (it->var.type & sbio::DELETION) {
							dlen = vit_->var.pos[0].begin - it->var.pos[0].begin;
							ilen = (vit_->var.pos[0].begin + vit_->var.pos[0].end) - (it->var.pos[0].begin + it->var.pos[0].end);
							cvar.type = (0 < dlen ? sbio::DELETION : 0) | sbio::INSERTION;
							cvar.pos[0] = it->var.pos[0];
							cvar.pos[0].end = cvar.pos[0].begin + dlen - 1;
							cvar.pos[1] = it->var.pos[0];
							cvar.pos[1].begin = it->var.pos[0].end;
							cvar.pos[1].end = cvar.pos[1].begin + ilen - 1;
						}
						else {
							dlen = (vit_->var.pos[0].begin + vit_->var.pos[0].end) - (it->var.pos[0].begin + it->var.pos[0].end);
							ilen = vit_->var.pos[0].begin - it->var.pos[0].begin;
							cvar.type = (0 < dlen ? sbio::DELETION : 0) | sbio::INSERTION;
							cvar.pos[0] = it->var.pos[0];
							cvar.pos[0].begin = it->var.pos[0].end;
							cvar.pos[0].end = cvar.pos[0].begin + dlen - 1;
							cvar.pos[1] = it->var.pos[0];
							cvar.pos[1].end = cvar.pos[1].begin + ilen - 1;
						}
						cvar.alt = it->var.alt + "/" + vit_->var.alt;
						it->var = cvar;
						it->idx2 = vit_->idx2;
						vit_->var.type = 0;
					}
				}
			}
			else if (it->var.type & sbio::TRANSLOCATION && vit_->var.type & sbio::TRANSLOCATION &&
				it->var.pos[0].idx == vit_->var.pos[1].idx && it->var.pos[1].idx == vit_->var.pos[0].idx &&
				it->var.pos[0].dir == vit_->var.pos[1].dir && it->var.pos[1].dir == vit_->var.pos[0].dir) {
				if (it->var.type & sbio::INVERSION) {
					if (it->var.pos[0].dir) {
						if (vit_->var.pos[1].begin < it->var.pos[0].begin &&
							it->var.pos[1].begin <= vit_->var.pos[0].begin) {
							svar_data cvar;
							int dlen = it->var.pos[0].begin - vit_->var.pos[1].begin - 1,
								ilen = vit_->var.pos[0].begin - it->var.pos[1].begin + 1;
							cvar.type = (0 < dlen ? sbio::DELETION : 0) | sbio::TRANSLOCATION | sbio::INSERTION | sbio::INVERSION;
							cvar.pos[0] = vit_->var.pos[1];
							++cvar.pos[0].begin;
							cvar.pos[0].end = cvar.pos[0].begin + dlen - 1;
							cvar.pos[1] = it->var.pos[0];
							cvar.pos[1].end = cvar.pos[1].begin + ilen - 1;
							cvar.alt = vit_->var.alt + "/" + it->var.alt;
							it->var = cvar;
							it->idx2 = vit_->idx2;
							vit_->var.type = 0;
						}
					}
					else {
						if (it->var.pos[0].begin < vit_->var.pos[1].begin &&
							vit_->var.pos[0].begin <= it->var.pos[1].begin) {
							svar_data cvar;
							int dlen = vit_->var.pos[1].begin - it->var.pos[0].begin - 1,
								ilen = it->var.pos[1].begin - vit_->var.pos[0].begin + 1;
							cvar.type = (0 < dlen ? sbio::DELETION : 0) | sbio::TRANSLOCATION | sbio::INSERTION | sbio::INVERSION;
							cvar.pos[0] = it->var.pos[0];
							++cvar.pos[0].begin;
							cvar.pos[0].end = cvar.pos[0].begin + dlen - 1;
							cvar.pos[1] = vit_->var.pos[0];
							cvar.pos[1].end = cvar.pos[1].begin + ilen - 1;
							cvar.alt = it->var.alt + "/" + vit_->var.alt;
							it->var = cvar;
							it->idx2 = vit_->idx2;
							vit_->var.type = 0;
						}
					}
				}
				else {
					if (it->var.pos[0].dir) {
						if (vit_->var.pos[1].begin < it->var.pos[0].begin &&
							vit_->var.pos[0].begin <= it->var.pos[1].begin) {
							svar_data cvar;
							int dlen = it->var.pos[0].begin - vit_->var.pos[1].begin - 1,
								ilen = it->var.pos[1].begin - vit_->var.pos[0].begin + 1;
							cvar.type = (0 < dlen ? sbio::DELETION : 0) | sbio::TRANSLOCATION | sbio::INSERTION;
							cvar.pos[0] = vit_->var.pos[1];
							++cvar.pos[0].begin;
							cvar.pos[0].end = cvar.pos[0].begin + dlen - 1;
							cvar.pos[1] = vit_->var.pos[0];
							cvar.pos[1].end = cvar.pos[1].begin + ilen - 1;
							cvar.alt = vit_->var.alt + "/" + it->var.alt;
							it->var = cvar;
							it->idx2 = vit_->idx2;
							vit_->var.type = 0;
						}
					}
					else {
						if (it->var.pos[0].begin < vit_->var.pos[1].begin &&
							it->var.pos[1].begin <= vit_->var.pos[0].begin) {
							svar_data cvar;
							int dlen = vit_->var.pos[1].begin - it->var.pos[0].begin - 1,
								ilen = vit_->var.pos[0].begin - it->var.pos[1].begin + 1;
							cvar.type = (0 < dlen ? sbio::DELETION : 0) | sbio::TRANSLOCATION | sbio::INSERTION;
							cvar.pos[0] = it->var.pos[0];
							++cvar.pos[0].begin;
							cvar.pos[0].end = cvar.pos[0].begin + dlen - 1;
							cvar.pos[1] = it->var.pos[1];
							cvar.pos[1].end = cvar.pos[1].begin + ilen - 1;
							cvar.alt = it->var.alt + "/" + vit_->var.alt;
							it->var = cvar;
							it->idx2 = vit_->idx2;
							vit_->var.type = 0;
						}
					}
				}
			}
			else if (it->var.type == sbio::INVERSION && vit_->var.type == sbio::INVERSION &&
				it->var.pos[0].idx == vit_->var.pos[0].idx &&
				it->var.pos[0].dir == vit_->var.pos[1].dir) {
				if (it->var.pos[0].dir &&
					vit_->var.pos[1].begin < it->var.pos[0].begin && it->var.pos[1].begin <= vit_->var.pos[0].begin) {
					svar_data cvar;
					int dlen = it->var.pos[0].begin - vit_->var.pos[1].begin,
						ilen = vit_->var.pos[0].begin - it->var.pos[1].begin + 1;
					cvar.type = (0 < dlen ? sbio::DELETION : 0) | sbio::INVERSION | sbio::INSERTION;
					cvar.pos[0] = vit_->var.pos[1];
					++cvar.pos[0].begin;
					cvar.pos[0].end = cvar.pos[0].begin + dlen - 1;
					cvar.pos[1] = it->var.pos[1];
					cvar.pos[1].end = cvar.pos[1].begin + ilen - 1;
					cvar.alt = vit_->var.alt + "/" + it->var.alt;
					it->var = cvar;
					it->idx2 = vit_->idx2;
					vit_->var.type = 0;
				}
				else if (!it->var.pos[0].dir &&
					it->var.pos[0].begin < vit_->var.pos[1].begin && vit_->var.pos[0].begin <= it->var.pos[1].begin) {
					svar_data cvar;
					int dlen = vit_->var.pos[1].begin - it->var.pos[0].begin,
						ilen = it->var.pos[1].begin - vit_->var.pos[0].begin + 1;
					cvar.type = (0 < dlen ? sbio::DELETION : 0) | sbio::INVERSION | sbio::INSERTION;
					cvar.pos[0] = it->var.pos[0];
					++cvar.pos[0].begin;
					cvar.pos[0].end = cvar.pos[0].begin + dlen - 1;
					cvar.pos[1] = vit_->var.pos[0];
					cvar.pos[1].end = cvar.pos[1].begin + ilen - 1;
					cvar.alt = it->var.alt + "/" + vit_->var.alt;
					it->var = cvar;
					it->idx2 = vit_->idx2;
					vit_->var.type = 0;
				}
			}
		}*/

	}
}
inline void findVar(VAPair& vap, DNASeqTrie2& trie, Sequence& seq, SeqList& reference) {
	sfori(vap.aligns) {
		auto& align = vap.aligns[i];
		findSmallVar(i, align, seq, reference, vap.svariants);
		if (i < vap.aligns.size() - 1) findLargeVariant(i, i + 1, align, vap.aligns[i + 1], trie.queries[0], vap.lvariants);
	}
	if (!vap.lvariants.empty()) findComplexVariant(vap.lvariants);
}
inline void filterVar(VAPair& vap, float minQual = 50.f, Array<sregion>* target = nullptr, AnnotDB* adb = nullptr) {
	slib::sbio::VarParam varp;
	sforin(i, 0, 3) varp.vcp.min_qual[i] = minQual;
	auto count = vap.svariants.size();
	sfor(vap.svariants) {
		auto& variant = $_.second;
		if (variant.qual < minQual ||
			(target && !target->at(variant.pos[0].idx).overlap(variant.pos[0]))) {
			variant.type = 0; --count;
		}
		else {
			bool identified = false;
			if (adb) {
				identified = adb->verifyVariant(variant, varp);
			}
			if (!identified) variant.varid = "Unknown variant";
		}
	}
	vap.svariants.sort([](const Pair<int, Variant>& p1, const Pair<int, Variant>& p2) {
		if (p1.second.type == 0) return false;
		if (p2.second.type == 0) return true;
		return p1.second < p2.second;
		});
	vap.svariants.resize(count);


	/*
	if (!vap.lvariants.empty()) {
		size_t size = vr.lvariants.size();
		sforeach(vr.lvariants) {
			if (!E_.var.type) { --size; continue; }
			if ((E_.var.pos[0].length() < par->min_size && E_.var.pos[1].length() < par->min_size) ||
				!par->isTargeted(&it->var)) {
				E_.var.type = 0;
				--size; continue;
			}
		}
		std::sort(vr.lvariants.begin(), vr.lvariants.end(),
			[](const va_pair& vp1, const va_pair& vp2) {
				if (!vp1.var.type) return false;
				if (!vp2.var.type) return true;
				return vp1.idx1 < vp2.idx1;
			});
		vap.lvariants.resize(size);
	}
	*/
	//if (adb) adb->fetchVariant();
}
inline int convABI2SLIB(int s) {
	switch (s) {
	case 0:
		return 4;
	case 1:
		return 1;
	case 2:
		return 8;
	case 3:
		return 2;
	}
	return -1;
}
inline void baserecall(sbio::Sequence& query) {
	Pair<int, int> max_base;
	auto base_pos = reinterpret_cast<sshort*>(&query.attribute["PLOC"][1]["data"].data()[0]);
	sshort* base_raw[4] = {
		reinterpret_cast<sshort*>(&query.attribute["DATA"][8]["data"].data()[0]),
		reinterpret_cast<sshort*>(&query.attribute["DATA"][9]["data"].data()[0]),
		reinterpret_cast<sshort*>(&query.attribute["DATA"][10]["data"].data()[0]),
		reinterpret_cast<sshort*>(&query.attribute["DATA"][11]["data"].data()[0])
	};
	sforin(b, 0, query.length()) {
		auto bp = (int)base_pos[b];
		if (query[b] == 0x01 || query[b] == 0x02 || query[b] == 0x04 || query[b] == 0x08) {
			max_base = Pair<int, int>(0, base_raw[0][bp]);
			sforin(d, 1, 4) { if (max_base.second < base_raw[d][bp]) max_base = Pair<int, int>(d, base_raw[d][bp]); }
			sforin(d, 0, 4) {
				if (d == max_base.first) continue;
				if ((max_base.second / 2) < (int)base_raw[d][bp]) {
					query[b] |= (subyte)convABI2SLIB(d);
				}
			}
		}
	}
}
void _makeAlignStr(Sequence& query, SeqList& ref, AlignPair& al, String* buf) {
	buf[0] = al.alref(ref.raw(al.ref));
	buf[1] = al.match();
	buf[2] = al.alque(query.raw().substring(al.query.begin, al.query.length(true)));
}
void _showAlign(Sequence& query, SeqList& ref, AlignPair& al, String* buf, IOStream& ostream) {
	_makeAlignStr(query, ref, al, buf);
	if (buf[0].size() < FASTA_ROW_COUNT) {
		ostream << sstr::lfill(S(al.ref.dir ? al.ref.end + 1 : al.ref.begin + 1), ' ', 12) << " : " << buf[0] << LF <<
			SP * 15 << buf[1] << LF <<
			sstr::lfill(S(al.query.begin + 1), ' ', 12) << " : " << buf[2] << LF;
	}
	else {
		auto count = (buf[0].size() - 1) / FASTA_ROW_COUNT;
		auto label = srange(al.ref.dir ? al.ref.end + 1 : al.ref.begin + 1, al.query.begin + 1);
		sforin(i, 0, count) {
			sforin(s, 0, 3) buf[s + 3] = buf[s].substring(FASTA_ROW_COUNT * i, FASTA_ROW_COUNT);
			ostream << sstr::lfill(S(label.begin), ' ', 12) << " : " << buf[3] << LF <<
				SP * 15 << buf[4] << LF <<
				sstr::lfill(S(label.end), ' ', 12) << " : " << buf[5] << LF * 2;

			label.begin += (al.ref.dir ? -1 : 1) * (buf[3].size() - buf[3].count("-"));
			label.end += (buf[5].size() - buf[5].count("-"));
		}
		sforin(s, 0, 3) buf[s + 3] = buf[s].substring(FASTA_ROW_COUNT * count);
		ostream << sstr::lfill(S(label.begin), ' ', 12) << " : " << buf[3] << LF <<
			SP * 15 << buf[4] << LF <<
			sstr::lfill(S(label.end), ' ', 12) << " : " << buf[5] << LF;

	}
	ostream.flush();
}

inline void printVarAlign(VAPair& vap, Sequence& query, SeqList& ref, IOStream& ostream, bool show) {
	String buf[6];
	if (vap.svariants.empty() && vap.lvariants.empty()) {
		ostream.print("Variant : Not found.", LF);
		if (show) {
			sforeach(align, vap.aligns) {
				_showAlign(query, ref, align, buf, ostream);
				ostream << NL; ostream.flush();
			}
		}
	}
	else {
		//
		if (vap.svariants.size()) {
			vap.svariants.sort([](const Pair<int, Variant>& p1, const Pair<int, Variant>& p2) {
				return p1.first < p2.first;
				});
			//
			sfori(vap.aligns) {
				auto count = 0;
				sforeach(vpair, vap.svariants) {
					if (i < vpair.first) break;
					else if (vpair.first < i) continue;
					//
					auto& variant = vpair.second;
					ostream << "Variant : " << variant.varid << TAB <<
						ref[variant.pos[0].idx].name << ":" << variant.pos[0].begin << "-" <<
						(variant.type == sbio::INSERTION ? variant.alt + "-" : "") <<
						variant.pos[0].end << TAB << sbio::sutil::varTypeDesc(variant.type) << TAB;
					if (variant.type == sbio::SNV || variant.type == sbio::MNV) {
						ostream << ref[variant.pos[0].idx].raw(variant.pos[0].begin - 1, variant.pos[0].length(true)) << "/" << variant.alt << TAB;
					}
					if (variant.genotype == HETERO_VAR) ostream << "(Hetrozygous)";
					++count;
				}
				if (show && 0 < count) {
					_showAlign(query, ref, vap.aligns[i], buf, ostream);
					ostream << NL; ostream.flush();
				}
				ostream << NL; ostream.flush();
			}
		}
		//
		if (vap.lvariants.size()) {
			sforeach(vpair, vap.lvariants) {
				auto& variant = vpair.second;
				ostream << "Variant : " << variant.varid << TAB;
				if (variant.type & INSERTION) {
					if (variant.pos[1].idx == -1) {
						ostream << ref[variant.pos[0].idx].name << ":" << variant.pos[0].begin << "-" <<
							variant.alt << "-" <<
							variant.pos[0].begin + 1 << TAB << sbio::sutil::varTypeDesc(variant.type);
					}
					else {
						ostream << ref[variant.pos[0].idx].name << ":" << variant.pos[0].begin << "- <" <<
							variant.alt.split("|")[0] <<
							ref[variant.pos[1].idx].name << ":" << variant.pos[1].begin << "-" <<
							variant.pos[1].end << "(" << (variant.pos[1].dir ? "-" : "+") << ")" << "> -" <<
							variant.alt.split("|")[1] <<
							variant.pos[0].end << TAB << sbio::sutil::varTypeDesc(variant.type) << LF * 2;


					}
				}
				else if (variant.pos[1].idx == -1) {
					ostream << ref[variant.pos[0].idx].name << ":" << variant.pos[0].begin << "-" <<
						(variant.alt.size() ? variant.alt + "-" : "") <<
						variant.pos[0].end << TAB << sbio::sutil::varTypeDesc(variant.type);
				}
				else {
					if (variant.type & TRANSLOCATION) {
						ostream << ref[variant.pos[0].idx].name << ":" << variant.pos[0].begin << "-" <<
							(variant.alt.size() ? variant.alt + "-" : "") <<
							ref[variant.pos[1].idx].name << ":" << variant.pos[1].begin <<
							"(" << (variant.pos[1].dir ? "-" : "+") << ")" << TAB <<
							sbio::sutil::varTypeDesc(variant.type);
					}
					else if (variant.type & INVERSION) {
						ostream << ref[variant.pos[0].idx].name << ":" << variant.pos[0].begin << "-" <<
							(variant.alt.size() ? variant.alt + "-" : "") <<
							ref[variant.pos[1].idx].name << ":" << variant.pos[1].begin <<
							"(" << (variant.pos[1].dir ? "-" : "+") << ")" << TAB <<
							sbio::sutil::varTypeDesc(variant.type);
					}
				}


				ostream << NL; ostream.flush();
			}
		}
	}
}
/*
inline void showVSResult(int res, Array<VAPair>& results, Sequence& query, SeqList& reference, bool show, IOStream& ostream) {
	switch (res) {
	case SHORT_SEQUENCE_ERROR:
	{
		ostream.print("!! Sequence is too short. !!");
		break;
	}
	case TOO_MUCH_N_ERROR:
	{
		ostream.print("!! There are too much N to analyze. !!");
		break;
	}
	case NOT_ALIGNED_ERROR:
	{
		ostream.print("!! Sequence is not aligned to the reference. !!");
		break;
	}
	case MULTI_HIT_ERROR:
	{
		ostream.print("* There are multiple candidates.");
		sforeach(result, results) {
			printVarAlign(result, query, reference, ostream, show);
			if (&result == &results[-1]) break;
			ostream.print(S("~") * 80);
		}
		break;
	}
	default:
		printVarAlign(results[0], query, reference, ostream, show);
		break;
	}
}
*/
Response& Moirei::varSearch(const SDictionary& pref) {
	SeqList reference;
	reference.load(pref["reference"]);
	//
	Array<sregion> target;
	if (pref.hasKey("target")) {
		slib::sbio::BEDFile bed(pref["target"]);
		bed.setRef(reference);
		bed >> target;
		sfor(target) $_.shift(-1);
	}
	//
	if (param.oformat == "auto") param.oformat = "txt";
	bool jsonout = param.oformat = "json";
	bool showalign = pref["alignment"];
	bool showchrom = pref["chromatogram"];
	//
	AnnotDB adb;
	bool annot = false;
	if (pref.hasKey("annotdb")) {
		adb.open(pref["annotdb"]);
		annot = true;
	}
	//
	SeqSearchParam seqp(DNA_SEQ4);
	//if (pref.hasKey("gap")) seqp.max_gap = pref["gap"]; else seqp.max_gap = 2;
	//if (pref.hasKey("miss")) seqp.max_miss = pref["miss"]; else seqp.max_miss = 2;
	seqp.max_gap = 0;
	seqp.max_miss = 0;
	if (pref.hasKey("match")) seqp.min_match = pref["match"]; else seqp.min_match = 16;
	if (pref.hasKey("seed")) seqp.setSeed(pref["seed"]);
	if (pref.hasKey("threshold")) seqp.ext_threshold = pref["threshold"]; else seqp.ext_threshold = 0.75;
	//
	DNASeqTrie2 trie(&seqp);
	SeqSearch search(&seqp);
	Array<AlignPair*> aligns;
	Array<Array<Pair<srange, salign*>>> maps;
	Array<VAPair> results;
	//
	float min_base_qual = 50.f;
	if (pref.hasKey("min-baseq")) min_base_qual = pref["min-baseq"];
	//
	bool recall = pref["recall"];

	//
	sfor(pref["_args_"]) {
		// Input
		auto input = sfs::absolutePath($_);
		if (sfs::exist(input)) {
			// Output
			if (param.outdir.empty()) param.ostream.setStdOStream(std::cout);
			else {
				param.ofile.open(sfs::joinPath(param.outdir, sfs::fileName(input, false) + ".txt"), MAKE);
				param.ostream.setFileOStream(param.ofile);
			}
		}
		else {
			SPrint("Input '", $_, "' is not found.");
			//return NOT_FOUND_ERROR;
		}
		Sequence query(DNA_SEQ);
		query.load(input);
		if (query.length() < seqp.seed) {
			//showVSResult(SHORT_SEQUENCE_ERROR, results, query, reference, showalign, param.ostream);
			//return 0;
		}
		// Re-basecall
		if (recall) baserecall(query);
		//
		trie.addQuery(query);
		trie.complete();
		search.search(reference, trie);
		//
		aligns.clear();
		sfor(search.aligns) {
			sforeach(align, $_) aligns.add(&align);
		}
		//
		if (aligns.empty()) {
			int count = 0;
			sforeach(base, query) { if (base == 0 || base == 15) ++count; }
			//if (query.length() / 5 < count) showVSResult(TOO_MUCH_N_ERROR, results, query, reference, showalign, param.ostream);
			//else showVSResult(NOT_ALIGNED_ERROR, results, query, reference, showalign, param.ostream);
			//return 0;
		}
		//
		aligns.sort([](const AlignPair* a1, const AlignPair* a2) { return a2->score < a1->score; });
		auto max_score = aligns[0]->score;
		sfor(aligns) {
			if ($_->score < max_score) break;
			Array<Pair<srange, salign*>> map;
			append2map($_, map, query.length());
			fillmap(map, Range<ArrayIterator<AlignPair*>>($ + 1, aligns.end()), maps, query.length(), &seqp);
		}
		//
		results.resize(maps.size());
		sfor2(maps, results) map2res($_1, $_2);
		results.sort([](const VAPair& p1, const VAPair& p2) {
			return p1.score() > p2.score();
			});
		//
		max_score = results[0].score();
		auto count = 0;
		sfor(results) {
			if ($_.score() < max_score) break;
			findVar($_, trie, query, reference);
			filterVar($_, 20.f, (target.empty() ? nullptr : &target), (adb.isOpened() ? &adb : nullptr));
			++count;
		}
		//
		results.resize(count);
		//if (1 < count) showVSResult(MULTI_HIT_ERROR, results, query, reference, showalign, param.ostream);
		//else showVSResult(ANALYSIS_COMPLETED, results, query, reference, showalign, param.ostream);
		// Close filestream if need.
		if (param.ofile.isOpened()) param.ofile.close();
	}
	return param.response;
}
