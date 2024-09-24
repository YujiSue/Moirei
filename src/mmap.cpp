#include "moirei.h"
Response& Moirei::motifMap(const SDictionary& pref) {
	try {
		param.setPref(pref);
		//
		float X_OFFSET = 5.f,
			Y_OFFSET = 50.f,
			PROT_BAR_OFFSET = 25.f,
			LSCALE_HEIGHT = 10.f,
			BOX_HEIGHT = 5.f,
			VERTICAL_SHIFT = 7.5f,
			SECTOR_MARGIN = 30.f,
			LABEL_MARGIN = 2.5f,
			MAIN_TEXT_SIZE = 12.f,
			SUB_TEXT_SIZE = 10.f,
			NOTE_TEXT_SIZE = 8.f;

		//
		AnnotDB adb(pref["annotdb"]);

		//
		auto& refIndex = adb.chrindex;
		//
		TextAttribute boxAttr;
		boxAttr.size = SUB_TEXT_SIZE;
		//
		auto width = pref["width"] == "auto" ? 640 : pref["width"].intValue();
		SCanvas mmap(width, 1);
		//
		float ratio = 1.0f;
		//
		int pid;
		GeneInfo ginfo;
		TranscriptInfo tinfo;
		AnnotInfo pinfo;
		//
		String query = (const char*)pref["_args_"][0];
		if (pref.hasKey("query-type")) {
			if (pref["query-type"] == "gene") {
				auto& idx = adb.searchGenes(query, MATCH::EXACT, {
					D_("search-synonym", true),
					D_("search-xref", true)
					});
				if (idx.empty()) throw NotFoundException(nofoundErrorText(query, "Gene database"));
				ginfo = adb.geneInfo(idx[0], { D_("transcript", true) });
				auto &proteins = adb.proteinsOf(idx[0]);
				if (proteins.empty()) throw NotFoundException(nofoundErrorText(query, "Protein database"));
				pinfo = proteins[0];
				sforeach(rna, ginfo.transcripts) {
					

				}

			}
			else if (pref["query-type"] == "transcript") {
				auto& transcripts = adb.getTranscripts(query);
				/*
				if (transcripts.size()) pos = genes[0];
				else throw NotFoundException(nofoundErrorText(query, "gene name/symbol/ID in annotation database."));
				pos = genes[0];
				*/
			}
		}
		else {
			/*
			auto& proteins = adb.getProteins(query);
			if (genes.size()) pos = genes[0];
			else throw NotFoundException(nofoundErrorText(query, "gene name/symbol/ID in annotation database."));
			pos = genes[0];
			int margin = pos.length();
			margin /= 10;
			if (10 < margin) margin -= (margin % 10);
			pos.begin -= margin;
			pos.end += margin;
			*/
		}
		ratio = (float)(width - X_OFFSET * 2) / (float)pinfo.length(true);
		//
		sushort display = (sushort)ANNOT_CATEGORY::MOTIF;
		if (pref.hasKey("mapping")) {
			display |= (sushort)ANNOT_CATEGORY::USER_DEFINED;
		}		
		// Protein
		{
			SFigure protArea;
			protArea.setId("Protein");
			SFigure protLine = SLine(svec2f(X_OFFSET, PROT_BAR_OFFSET), svec2f(width- X_OFFSET, PROT_BAR_OFFSET), Stroke(sstyle::SOLID_LINE, 2.f)),
				plotLabel = SCaption(X_OFFSET, MAIN_TEXT_SIZE, pinfo.name);

			SFigure scaleLabels;
			Path2D<float> longScale;
			//
			int lbin = smath::power(10, log(pinfo.length()) / log(10.0));
			if (pinfo.length() / lbin < 3) lbin /= 2;
			//
			int lsc = (pinfo.begin / lbin + 1) * lbin;
			//
			TextAttribute sclAttr;
			sclAttr.size = NOTE_TEXT_SIZE;
			//
			while (lsc < pinfo.end) {
				auto scLabel = SCaption((float)(lsc - pinfo.begin) * ratio - sclAttr.size / 2.f, PROT_BAR_OFFSET - LSCALE_HEIGHT / 2.f - NOTE_TEXT_SIZE / 2.f, S(lsc), sclAttr);
				scaleLabels.draw(scLabel);
				longScale.moveTo(svec2f((float)(lsc - pinfo.begin) * ratio, PROT_BAR_OFFSET - LSCALE_HEIGHT / 2.f));
				longScale.lineTo(svec2f((float)(lsc - pinfo.begin) * ratio, PROT_BAR_OFFSET + LSCALE_HEIGHT / 2.f));
				lsc += lbin;
			}
			protArea
				.draw(protLine)
				.draw(plotLabel)
				.draw(scaleLabels)
				.draw(SPath(longScale));
			mmap.draw(protArea);
		}
		if (display & (sushort)ANNOT_CATEGORY::MOTIF) {
			SFigure motifArea;
			sobj sopt;
			if (pref.hasKey("motif-source")) sopt = { D_("source", pref["motif-source"]) };
			auto& motifs = adb.motifsOf(pinfo.record, sopt);
			sforeach(motif, motifs) {
				SFigure motfig,
					motlabel = SCaption((motif.begin - pinfo.begin) * ratio + X_OFFSET, Y_OFFSET, motif.name, boxAttr),
					motbox = SRectangle((motif.begin - pinfo.begin) * ratio + X_OFFSET, Y_OFFSET + LABEL_MARGIN, (float)motif.length(true) * ratio, BOX_HEIGHT);
				motbox.setFillColor(scolor::BLUE);
				motfig
					.draw(motlabel)
					.draw(motbox);
				motfig.setId(motif.name);
				motifArea.draw(motfig);
			}
			if (motifs.size()) {
				moir::rearrange(motifArea, VERTICAL_SHIFT);
				mmap.draw(motifArea);
				Y_OFFSET = mmap.boundary().ori_y + mmap.boundary().height + SECTOR_MARGIN;
			}
		}

		if (display & (sushort)ANNOT_CATEGORY::USER_DEFINED) {
			SFigure customArea;
			auto objs = pref["mapping"].split(",");
			sforeach(args, objs) {
				auto obj = args.parse(";", "=");
				SFigure customFig;
				auto start = obj["pos"].intValue(), len = obj["length"].intValue();
				srange range = srange(start, start + len - 1);
				if (!obj.hasKey("type")) obj["type"] = "pos";
				if (obj["type"] == "genome") {
					range = sbio::sutil::codingSite(range, &tinfo, tinfo.gene->dir);
					range.begin /= 3; range.end /= 3;
				}
				else if (obj["type"] == "transcript") {
					range.begin /= 3; range.end /= 3;
				}
				Color col = obj.hasKey("color") ? Color(obj["color"]) : scolor::RED;
				SFigure customLabel = SCaption((range.begin - pinfo.begin) * ratio + X_OFFSET, Y_OFFSET, obj["label"]),
					customRect = SRectangle((range.begin - pinfo.begin) * ratio, Y_OFFSET + LABEL_MARGIN, range.length(true) * ratio, BOX_HEIGHT);
				customRect.setFillColor(col);
				customFig.draw(customLabel).draw(customRect);
				customArea.draw(customFig);
				moir::rearrange(customArea, VERTICAL_SHIFT);
			}
			mmap.draw(customArea);
		}
		mmap.resize(width, mmap.boundary().height + SECTOR_MARGIN);
		if (param.oformat == "_obj_") param.response.attribute["result"] = mmap;
		if (param.output.empty()) {
			param.ostream.setStrOStream(param.response.output);
			param.ostream.print(sxml::svgNode(mmap)->toString());
		}
		else mmap.save(param.output);
		if (param.response.output) SPrint(param.response.output);
	}
	catch (Exception ex) {
		ex.print();
		param.response = ex;
	}
	return param.response;

}



//extern inline void makeMotifMap(SBSeqList &reference, SBAnnotDB &annotdb, const sbpos &range, Array<>) {
extern inline void getMotifs() {

}

//}
/*
int main() {
	slib::String dir = "/content/protein/worm/fasta", seq, rseq, aseq;
	ubytearray rseqi, aseqi;
	slib::sbio::SBSeqList ref;
	ref.load("/content/nematode.bin");
	slib::sbio::SBAnnotDB annot;
	annot.open("/content/nematode.db");
	annot.load(slib::sbio::SBAnnotDB::LOAD_GENE | slib::sbio::SBAnnotDB::LOAD_TRANS);
	sforeach_(git, annot.genes) {
		if (git->type != slib::sbio::PROTEIN_CODING) continue;
		sforeach_(tit, git->transcripts) {
			if ((*tit)->type != slib::sbio::M_RNA) continue;
			seq.clear();
			sforeach_(sit, (*tit)->structures) {
				if (sit->type == slib::sbio::CDS) {
					seq += ref[sit->idx]->raw(sit->begin - 1, sit->length(true));
				}
			}
			if (git->dir) seq = slib::sbio::sseq::dcompseq(seq);
			rseq.resize(seq.size(), '\0');
			slib::sbio::sseq::dtrans((const subyte*)seq.cstr(), 0, seq.size(), (subyte*)&rseq[0]);
			rseqi.resize(rseq.size());
			slib::sbio::sseq::rencode((const subyte*)rseq.cstr(), 0, rseq.size(), &rseqi[0]);
			aseqi.resize(rseqi.size() / 3, 0);
			slib::sbio::sseq::atrans(rseqi.ptr(), 0, rseqi.size(), &aseqi[0], slib::sbio::DEFAULT_CODON);
			aseq.resize(aseqi.size() - 1, '\0');
			slib::sbio::sseq::adecode(aseqi.ptr(), 0, aseqi.size() - 1, (subyte*)&aseq[0]);
			sio::SFile file(dir + PATH_SEPARATOR + (*tit)->name + ".fa", sio::CREATE);
			file << ">" << (*tit)->name << LF;
			file << aseq;
			file.flush();
			file.close();
		}
	}

	Map<String, Array<srange>> targets, converted;
	targets["spv-1"] = Array<srange>({ srange(7915067, 7915643), srange(7925596, 7925614) });

	slib::sbio::SBSeqList ref;
	ref.load("/content/nematode.bin");
	slib::sbio::SBAnnotDB annot;
	annot.open("/content/nematode.db");
	annot.load(slib::sbio::SBAnnotDB::LOAD_GENE | slib::sbio::SBAnnotDB::LOAD_TRANS);
	sbio::SBAnnotDB::geneparray genes;
	sforeach_(git, targets) {
		annot.geneInfo(genes, git->key);
		auto& region = git->value;
		sforeach_(tit, genes[0]->transcripts) {
			if ((*tit)->type != slib::sbio::M_RNA) continue;
			size_t len = 0;
			srange tpos[2];
			tpos[0].begin = -1; tpos[0].end = -1;
			tpos[1].begin = -1; tpos[1].end = -1;
			sforeach_(sit, (*tit)->structures) {
				if (sit->type == slib::sbio::CDS) {
					if (tpos[0].begin == -1) {
						if (region[0].begin <= sit->begin) tpos[0].begin = len;
						else if (region[0].begin <= sit->end) tpos[0].begin = len + region[0].begin - sit->begin;
					}
					if (tpos[0].end == -1) {
						if (region[0].end < sit->begin) tpos[0].end = len - 1;
						else if (region[0].end <= sit->end) tpos[0].begin = len + region[0].end - sit->begin;
					}
					if (tpos[1].begin == -1) {
						if (region[1].begin <= sit->begin) tpos[1].begin = len;
						else if (region[1].begin <= sit->end) tpos[1].begin = len + region[1].begin - sit->begin;
					}
					if (tpos[1].end == -1) {
						if (region[1].end < sit->begin) tpos[1].end = len - 1;
						else if (region[1].end <= sit->end) tpos[1].begin = len + region[1].end - sit->begin;
					}
					len += sit->length(true);
				}
			}
			if (genes[0]->dir) {
				tpos[0].begin = len - tpos[0].begin - 1;
				tpos[0].end = len - tpos[0].end - 1;
				tpos[1].begin = len - tpos[1].begin - 1;
				tpos[1].end = len - tpos[1].end - 1;
			}
			converted[(*tit)->name] = Array<srange>({ tpos[0], tpos[1] });
			String output;
			SSystem::exec("cat /content/drive/MyDrive/interpro.tsv | grep " + (*tit)->name + " | grep Pfam", output);
			auto list = output.split("\t");
			SPrint(list);

		}

	}

	return 0;
}

void saveSVG(const char* gene, const char* name, size_t len, stringarray domains, Map<String, srange>&& vars) {
	smedia::SCanvas canvas;
	float offsetx = 10.f, offsety = 25.f,
		width = 640.f, scale = (width - 2 * offsetx) / len;

	sfig aaseq(sshape::GROUP);
	aaseq->attribute()["name"] = name;

	sline aa_line(v2f(offsetx, offsety), v2f(width - offsetx, offsety));
	aa_line->setStrokeWidth(2.0f);
	//scalligraphy aa_label(offsetx, 12.0, name);
	//aa_label->setFont("Arial", 12.0f);
	aaseq->addFigure(aa_line);
	canvas.drawFigure(aaseq);
	//aaseq->addFigure(aa_label);
	//sfig aa_scale(sshape::GROUP);
	//smedia::SPath2D scales;
	sfig domain_area(sshape::GROUP);
	domain_area->attribute()["name"] = "domains";
	sforeach(domains) {
		auto values = E_.split("\t");
		//SPrint(values);
		//0:Transcript name, 1:hash?, 2:# of residues, 3:Method=Pfam, 4:PfamID?, 5:Domain name  6:start, 7:end, 8:qual, 9:?, 10:date?, 11:?, 12:Description
		sfig domain_g(sshape::GROUP);
		domain_g->attribute()["name"] = values[5];
		srect domain_rect(offsetx + scale * values[6].integer(), offsety - 7.5f, scale * (values[7].integer() - values[6].integer() + 1), 15.0f);
		domain_rect->setFillColor(smedia::ColorMap["blue"]);
		scalligraphy domain_label(domain_rect->boundary().ori_x, 16.0f, values[5]);
		domain_g->addFigure(domain_rect);
		domain_g->addFigure(domain_label);
		domain_area->addFigure(domain_g);
	}
	canvas.drawFigure(domain_area);
	offsety += 37.5f;
	sfig variant_area(sshape::GROUP);
	variant_area->attribute()["name"] = "variants";
	sforeach(vars) {
		sfig var_g(sshape::GROUP);
		var_g->attribute()["name"] = E_.key;
		srect var_rect(offsetx + scale * E_.value.begin, offsety, scale * E_.value.length(true), 10.f);
		var_rect->setFillColor(smedia::ColorMap["red"]);
		scalligraphy var_label(var_rect->boundary().ori_x, offsety + 25.f, E_.key);
		var_g->addFigure(var_rect);
		var_g->addFigure(var_label);
		variant_area->addFigure(var_g);
	}
	canvas.drawFigure(variant_area);
	canvas.resize(canvas.root().boundary().ori_x + canvas.root().boundary().width + 5.0f, canvas.root().boundary().ori_y + canvas.root().boundary().height + 5.0f);
	String dir = "/content/domain";
	SSystem::exec("mkdir -p " + dir);
	dir << "/" << gene;
	SSystem::exec("mkdir -p " + dir);
	dir << "/" << name << ".svg";
	canvas.save(dir);
}

int main2() {
	try {
		Map<String, Map<String, srange>> targets;
		targets["spv-1"] = Map<String, srange>();
		targets["spv-1"]["ok1498"] = srange(7915067, 7915643);
		targets["spv-1"]["tm16697"] = srange(7925596, 7925614);

		targets["dagl-2"] = Map<String, srange>();
		targets["dagl-2"]["tm2908"] = srange(782619, 782823);
		targets["dagl-2"]["tm16721"] = srange(778921, 778932);

		targets["bcc-1"] = Map<String, srange>();
		targets["bcc-1"]["tm3821"] = srange(11088640, 11089158);
		targets["bcc-1"]["tm10343"] = srange(11087679, 11087722);

		targets["dgk-3"] = Map<String, srange>();
		targets["dgk-3"]["gk110"] = srange(9165061, 9165329);
		targets["dgk-3"]["tm16695"] = srange(9166791, 9166820);

		targets["gcy-9"] = Map<String, srange>();
		targets["gcy-9"]["tm2816"] = srange(11333401, 11333784);
		targets["gcy-9"]["tm7632"] = srange(11334012, 11334057);

		targets["unc-73"] = Map<String, srange>();
		targets["unc-73"]["ev802"] = srange(4002535, 4004506);
		targets["unc-73"]["tm7912"] = srange(3999596, 3999633);

		targets["frk-1"] = Map<String, srange>();
		targets["frk-1"]["ok760"] = srange(10043707, 10045179);
		targets["frk-1"]["tm10434"] = srange(10045057, 10045101);

		targets["lad-2"] = Map<String, srange>();
		targets["lad-2"]["hd31"] = srange(2965250, 2965498);
		targets["lad-2"]["tm16720"] = srange(2973788, 2973800);

		slib::sbio::SBSeqList ref;
		ref.load("/content/nematode.bin");
		slib::sbio::SBAnnotDB annot;
		annot.open("/content/nematode.db");
		annot.load(slib::sbio::SBAnnotDB::LOAD_GENE | slib::sbio::SBAnnotDB::LOAD_TRANS);
		sbio::SBAnnotDB::geneparray genes;
		sforeach_(git, targets) {
			annot.geneInfo(genes, git->key);
			auto& vmap = git->value;
			auto vkeys = vmap.keyset();
			srange region[2] = { vmap[vkeys[0]], vmap[vkeys[1]] };
			sforeach_(tit, genes[0]->transcripts) {
				if ((*tit)->type != slib::sbio::M_RNA) continue;
				size_t len = 0;
				srange tpos[2];
				tpos[0].begin = -1; tpos[0].end = -1;
				tpos[1].begin = -1; tpos[1].end = -1;
				sforeach_(sit, (*tit)->structures) {
					if (sit->type == slib::sbio::CDS) {
						if (tpos[0].begin == -1) {
							if (region[0].begin <= sit->begin) tpos[0].begin = len;
							else if (region[0].begin <= sit->end) tpos[0].begin = len + region[0].begin - sit->begin;
						}
						if (tpos[0].end == -1) {
							if (region[0].end < sit->begin) tpos[0].end = len - 1;
							else if (region[0].end <= sit->end) tpos[0].end = len + region[0].end - sit->begin;
						}
						if (tpos[1].begin == -1) {
							if (region[1].begin <= sit->begin) tpos[1].begin = len;
							else if (region[1].begin <= sit->end) tpos[1].begin = len + region[1].begin - sit->begin;
						}
						if (tpos[1].end == -1) {
							if (region[1].end < sit->begin) tpos[1].end = len - 1;
							else if (region[1].end <= sit->end) tpos[1].end = len + region[1].end - sit->begin;
						}
						SPrint(tpos[0].begin, ",", tpos[0].end, " : ", tpos[1].begin, ",", tpos[1].end);
						len += sit->length(true);
					}
				}
				if (genes[0]->dir) {
					auto tmp1 = len - tpos[0].end - 1;
					auto tmp2 = len - tpos[0].begin - 1;
					tpos[0].begin = tmp1;
					tpos[0].end = tmp2;

					tmp1 = len - tpos[1].end - 1;
					tmp2 = len - tpos[1].begin - 1;
					tpos[1].begin = tmp1;
					tpos[1].end = tmp2;
				}
				SPrint(tpos[0].begin, ",", tpos[0].end, " : ", tpos[1].begin, ",", tpos[1].end);

				String output;
				SSystem::exec("cat /content/drive/MyDrive/interpro.tsv | grep " + (*tit)->name + " | grep -e Pfam -e SUPERFAMILY", output);
				auto lines = output.trimming().splitline();
				auto pro_len = len / 3 - 1;
				tpos[0].begin = tpos[0].begin / 3;
				tpos[0].end = tpos[0].end / 3;
				tpos[1].begin = tpos[1].begin / 3;
				tpos[1].end = tpos[1].end / 3;
				saveSVG(git->key, (*tit)->name, pro_len, lines, { kvpair<String, srange>(vkeys[0], tpos[0]), kvpair<String, srange>(vkeys[1], tpos[1]) });
			}
		}
	}
	catch (SException ex) { ex.print(); }
	return 0;
}
	*/
