#include "moirei.h"

inline void aaSubstitution(const SDictionary& pref) {
	Response res;
	AnnotDB adb;
	adb.open(pref["annotdb"]);




}

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
		range.end <map[0].first.begin || 
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
inline void makeVar(Cigar& cigar, salign* align, Sequence& seq, SeqList &reference, int &rpos, int &qpos, Variant& var) {
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

inline void printVarAlign(VAPair &vap, SeqList &ref, IOStream &ostream, const char *form, bool align) {
	String f(form);
	sforeach(vpair, vap.svariants) {
		auto& variant = vpair.second;
		if (f == "vcf") svcf::writeVariant(variant, ref, ostream);
		else {
			


		}
		if (align) {


		}
		ostream << NL; ostream.flush();
	}
	sforeach(vpair, vap.lvariants) {
		auto& variant = vpair.second;
		if (f == "vcf") svcf::writeVariant(variant, ref, ostream);
		else ostream << variant.toString(form) << TAB << variant.qual << NL; ostream.flush();
		if (align) {


		}
		ostream << NL; ostream.flush();
	}



	/*
	auto& v2 = res.lvariants;
	if (!v2.empty()) {
		sforeach_(vit, v2) {
			if (vit->var.type & sbio::INSERTION) {
				if (vit->var.type & sbio::TRANSLOCATION) {
					if (vit->var.type & sbio::INVERSION) {
						std::cout << par.ref[vit->var.pos[0].idx]->name << ":" << vit->var.pos[0].begin << "-<" <<
							par.ref[vit->var.pos[1].idx]->name << ":" << vit->var.pos[1].end << "-" << vit->var.pos[1].begin <<
							">-" << vit->var.pos[0].end << " (" << vit->var.pos[1].length() << " bp inverted translocational insertion)";
					}
					else {
						std::cout << par.ref[vit->var.pos[0].idx]->name << ":" << vit->var.pos[0].begin << "-<" <<
							par.ref[vit->var.pos[1].idx]->name << ":" << vit->var.pos[1].begin << "-" << vit->var.pos[1].end <<
							">-" << vit->var.pos[0].end << " (" << vit->var.pos[1].length() << " bp translocational insertion)";
					}
				}
				else if (vit->var.type & sbio::INVERSION) {
					std::cout << par.ref[vit->var.pos[0].idx]->name << ":" << vit->var.pos[0].begin << "-<" <<
						par.ref[vit->var.pos[1].idx]->name << ":" << vit->var.pos[1].end << "-" << vit->var.pos[1].begin <<
						">-" << vit->var.pos[0].end << " (" << vit->var.pos[1].length() << " bp inverted insertion)";
				}
				else if (vit->var.type & sbio::DELETION) {
					std::cout << par.ref[vit->var.pos[0].idx]->name << ":" << vit->var.pos[0].begin << "-<" <<
						par.ref[vit->var.pos[1].idx]->name << ":" << vit->var.pos[1].begin << "-" << vit->var.pos[1].end <<
						">-" << vit->var.pos[0].end << " (" << vit->var.pos[0].length() << " deletion + " << vit->var.pos[1].length() << " bp insertion)";
				}
				else {
					std::cout << par.ref[vit->var.pos[0].idx]->name << ":" << vit->var.pos[0].begin << "-<" <<
						par.ref[vit->var.pos[1].idx]->name << ":" << vit->var.pos[1].begin << "-" << vit->var.pos[1].end <<
						">-" << vit->var.pos[0].end << " (" << vit->var.pos[1].length() << " bp insertion)";
				}
			}
			else {
				if (vit->var.type & sbio::TRANSLOCATION) {
					if (vit->var.type & sbio::INVERSION) {
						std::cout << par.ref[vit->var.pos[0].idx]->name << ":" << vit->var.pos[0].begin << (vit->var.pos[0].dir ? "(-)" : "(+)") << "-" <<
							par.ref[vit->var.pos[1].idx]->name << ":" << vit->var.pos[1].begin << (vit->var.pos[1].dir ? "(-)" : "(+)") << " (translocational inversion)." << std::endl;
					}
					else {
						std::cout << par.ref[vit->var.pos[0].idx]->name << ":" << vit->var.pos[0].begin << "-" <<
							par.ref[vit->var.pos[1].idx]->name << ":" << vit->var.pos[1].begin << " (translocation)." << std::endl;
					}
				}
				else if (vit->var.type & sbio::INVERSION) {
					std::cout << par.ref[vit->var.pos[0].idx]->name << ":" << vit->var.pos[0].begin << (vit->var.pos[0].dir ? "(-)" : "(+)") << "-" <<
						par.ref[vit->var.pos[1].idx]->name << ":" << vit->var.pos[1].begin << (vit->var.pos[1].dir ? "(-)" : "(+)") << " (inversion)." << std::endl;
				}
				else if (vit->var.type & sbio::DELETION) {
					std::cout << par.ref[vit->var.pos[0].idx]->name << ":" << vit->var.pos[0].begin << "-" <<
						vit->var.pos[0].end << " (" << vit->var.pos[0].length() + 1 << " bp deletion)." << std::endl;
				}
				else if (vit->var.type & sbio::DUPLICATION) {
					std::cout << par.ref[vit->var.pos[0].idx]->name << ":" << vit->var.pos[0].begin << "-" <<
						vit->var.pos[0].end << " (" << vit->var.pos[0].length() + 1 << " bp duplication)." << std::endl;
				}
			}
		}
	}
	*/
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

inline void findSmallVar(int i, salign& al, Sequence &seq, SeqList& reference, Array<Pair<int, Variant>>& variants) {
	int rpos = al.ref.dir ? al.ref.end : al.ref.begin, qpos = al.query.begin;
	sfor(al.cigars) {
		if ($_.option == scigar::MMATCH ||
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
	if (al1.ref.dir) sbio::sutil::reverse(rng1, que.size());
	if (al2.ref.dir) sbio::sutil::reverse(rng2, que.size());
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
inline void findVar(VAPair& vap, DNASeqTrie2& trie, Sequence &seq, SeqList &reference) {
	sfori(vap.aligns) {
		auto& align = vap.aligns[i];
		findSmallVar(i, align, seq, reference, vap.svariants);
		if (i < vap.aligns.size() - 1) findLargeVariant(i, i + 1, align, vap.aligns[i + 1], trie.queries[0], vap.lvariants);
	}
	if (!vap.lvariants.empty()) findComplexVariant(vap.lvariants);
}
inline void filterVar(VAPair& vap, float minQual = 50.f, Array<sregion> *target = nullptr, AnnotDB *adb = nullptr) {
	/*
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
				auto& vars = adb->searchVariant(variant);
				if (vars.size()) {
					variant.varid = vars[0].name; identified = true;
				}
			}
			if(!identified) variant.varid = "Unknown variant";
		}
	}
	vap.svariants.sort([](const Pair<int, Variant> &p1, const Pair<int, Variant>& p2) {
		if (p1.second.type == 0) return false;
		if (p2.second.type == 0) return true;
		return p1.second < p2.second;
		});
	vap.svariants.resize(count);
	*/
	

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
inline subyte sangarAnalysis(const char *path, Array<VAPair> &results, SeqList &reference, Array<sregion> &target, AnnotDB &adb, SeqSearchParam &par) {
	if (!path || !sfs::exist(path)) return NO_SAMPLE_ERROR;
	Sequence query(DNA_SEQ);
	query.read(path);
	int len = query.length();
	if (len < par.seed) return SHORT_SEQUENCE_ERROR;
	//
	DNASeqTrie2 trie(&par);
	SeqSearch search(&par);
	trie.addQuery(query);
	trie.complete();
	search.search(reference, trie);
	//
	Array<AlignPair*> aligns;
	Array<Array<Pair<srange, salign*>>> maps;
	sfor(search.aligns) {
		sforeach(align, $_) aligns.add(&align);
	}
	//
	if (aligns.empty()) {
		int count = 0;
		sforeach(base, query) { if (base == 0 || base == 15) ++count; }
		if (len / 5 < count) return TOO_MUCH_N_ERROR;
		else return NOT_ALIGNED_ERROR;
	}
	//
	aligns.sort([](const AlignPair* a1, const AlignPair* a2) { return a1->score > a2->score; });
	auto max_score = aligns[0]->score;
	sfor(aligns) {
		if ($_->ref.idx < 0) continue;
		if ($_->score < max_score) break;
		Array<Pair<srange, salign*>> map;
		append2map($_, map, len);
		fillmap(map, Range<ArrayIterator<AlignPair*>>($ + 1, aligns.end()), maps, len, &par);
	}
	results.resize(maps.size());
	sfor2(maps, results) map2res($_1, $_2);
	results.sort([](const VAPair &p1, const VAPair& p2) {
		return p1.score() > p2.score();
		});
	max_score = results[0].score();
	auto count = 0;
	sfor(results) {
		if ($_.score() < max_score) break;
		findVar($_, trie, query, reference);
		filterVar($_, 50.f, (target.empty() ? nullptr : &target), (adb.isOpened() ? &adb : nullptr));
		++count;
	}
	results.resize(count);
	if (1 < count) return MULTI_HIT_ERROR;
	return 0xFF;
}
Response moir::varSearch(const SDictionary& pref) {
	Response res;
	IOStream ostream;
	//
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
	String oformat;
	bool showalign;
	if (pref.hasKey("oformat")) oformat = pref["oformat"];
	if (pref.hasKey("show-align")) showalign = pref["show-align"];

	AnnotDB adb;
	if (pref.hasKey("annotdb")) adb.open(pref["annotdb"]);
	//
	SeqSearchParam seqp(DNA_SEQ4);
	if (pref.hasKey("gap")) seqp.max_gap = pref["gap"]; else seqp.max_gap = 2;
	if (pref.hasKey("miss")) seqp.max_miss = pref["miss"]; else seqp.max_miss = 2;
	if (pref.hasKey("match")) seqp.min_match = pref["match"]; else seqp.min_match = 16;
	if (pref.hasKey("seed")) seqp.setSeed(pref["seed"]);
	if (pref.hasKey("threshold")) seqp.ext_threshold = pref["threshold"]; else seqp.ext_threshold = 0.75;
	//
	float min_base_qual = 50.f;
	if (pref.hasKey("min-baseq")) min_base_qual = pref["min-baseq"];
	//
	auto files = pref["_args_"];
	sfor(files) {
		SFile output;
		String outpath;
		if (pref.hasKey("outdir")) outpath = sfs::joinPath(pref["outdir"], sfs::fileName($_, false) + ".txt");
		else outpath = sfs::joinPath(sfs::splitPath($_).first, sfs::fileName($_, false) + ".txt");
		output.open(outpath, MAKE);
		ostream.setFileOStream(output);
		//
		Array<VAPair> results;
		auto ret = sangarAnalysis($_, results, reference, target, adb, seqp);
		switch (ret) {
		case NO_SAMPLE_ERROR:
		{
			ostream.print("Sequence file is not found.");
			res.code = FILE_NOT_EXIST_ERROR;
			res.error = "Sequence file is not found.";
			break;
		}
		case SHORT_SEQUENCE_ERROR:
		{
			ostream.print("Sequence is too short.");
			//res.code = FORMAT_ERROR
			res.error = "Sequence is too short.";
			break;
		}
		case TOO_MUCH_N_ERROR:
		{
			ostream.print("There are too much N to analyze.");
			//res.code = LOW_QUALITY_DATA_ERROR;
			res.error = "There are too much N to analyze.";
			break;
		}
		case NOT_ALIGNED_ERROR:
		{
			ostream.print("Sequence is not aligned to the reference.");
			//res.code = REFERENCE_MISMATCH_ERROR;
			res.error = "Sequence is not aligned to the reference.";
			break;
		}
		case MULTI_HIT_ERROR:
		{
			ostream.print("There are multiple candidates.");
			sforeach(result, results) {
				printVarAlign(result, reference, ostream, oformat, showalign);
				if (&result == &results[-1]) break;
				SPrint(S("*") * 60);
			}
			break;
		}
		default:
			printVarAlign(results[0], reference, ostream, oformat, showalign);
			break;
		}
		output.close();
	}
	return res;
}

typedef int getFilterInfo(slib::SDictionary &info);
typedef int customVarFilter(slib::sbio::Variant *var, sobj opts);
Response moir::varFilter(moir::Param& param, const SDictionary& pref) {
	Response res;
	//
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
	String oformat = pref.hasKey("oformat") ? (pref["oformat"] == "auto" ? "tsv" : (const char*)pref["oformat"]) : "tsv";
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
		auto output = sfs::joinPath(path.first, path.second.replace(sfs::extension($_), "filtered." + oformat));
		vlist.save(output, "cols=VarID,Chr1,Pos1,Len1,Ref,Alt/Ins,Type,Qual,Freq,Genotype,Gene,Site,Mutation,Substitution,Filter");
		SPrint("Save result to '", output, "'");
		if (plio) { plgio.exec($_, &vlist, piopts); }
		vlist.clearAll();
	}
	return res;
}