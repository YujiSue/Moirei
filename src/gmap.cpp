#include "moirei.h"

using CM_ = Pair<String, Color>;
Response& Moirei::geneMap(const SDictionary& pref) {
	try {
		param.setPref(pref);
		//
		float X_OFFSET = 5.f,
			Y_OFFSET = 50.f,
			CHR_BAR_OFFSET = 25.f,
			LSCALE_HEIGHT = 10.f,
			SSCALE_HEIGHT = 5.f,
			BOX_HEIGHT = 5.f,
			VERTICAL_SHIFT = 7.5f,
			SECTOR_MARGIN = 30.f,
			LABEL_MARGIN = 2.5f,
			MAIN_TEXT_SIZE = 12.f,
			SUB_TEXT_SIZE = 10.f,
			NOTE_TEXT_SIZE = 8.f;

		//
		RefPos pos;
		AnnotDB adb(pref["annotdb"]);
		auto& refIndex = adb.chrindex;
		//
		auto width = pref["width"] == "auto" ? 640 : pref["width"].intValue();
		SCanvas gmap(width, 1);
		//
		float ratio = 1.0f;
		//
		float arrowwidth = 5.f;
		//
		TextAttribute boxAttr;
		boxAttr.size = SUB_TEXT_SIZE;

		//
		String query = (const char*)pref["_args_"][0];
		String qtype = pref.hasKey("query-type") ? pref["query-type"].string() : "pos";
		if (qtype == "gene") {
			auto& genes = adb.getGenes(query);
			if (genes.size()) pos = genes[0];
			else throw NotFoundException(nofoundErrorText(query, "gene name/symbol/ID in annotation database."));
			pos = genes[0];
			int margin = pos.length();
			margin /= 10;
			if (10 < margin) margin -= (margin % 10);
			pos.begin -= margin;
			pos.end += margin;
		}
		else if (qtype == "pos") pos = RefPos::toPos(query, refIndex);
		ratio = (float)width / (float)pos.length(true);

		//
		sushort display = 0;
		if (pref["display"] == "auto") {
			display = (sushort)ANNOT_CATEGORY::GENE | (sushort)ANNOT_CATEGORY::TRANSCRIPT | (sushort)ANNOT_CATEGORY::MUTATION;
		}
		else {
			sforc(pref["display"].string()) {
				if ($_ == "G") display |= (sushort)ANNOT_CATEGORY::GENE;
				else if ($_ == "T") display |= (sushort)ANNOT_CATEGORY::TRANSCRIPT;
				else if ($_ == "V") display |= (sushort)ANNOT_CATEGORY::MUTATION;
				//else if ($_ == "F") display |= (sushort)ANNOT_CATEGORY::FEATURE;
			}
		}
		if (pref.hasKey("mapping") && pref["mapping"]) {
			display |= (sushort)ANNOT_CATEGORY::USER_DEFINED;
		}
		// Color
		Map<String, Color> colorMap = {
			CM_("gene", scolor::NAVY),
			CM_("pseudo", scolor::DIMGRAY),
			CM_("me", scolor::BROWN),
			CM_("CDS+", scolor::MAGENTA),
			CM_("CDS-", scolor::LIME),
			CM_("UTR", scolor::LIGHTGRAY),
			CM_("trna", scolor::DARKGREEN),
			CM_("rrna", scolor::DARKRED),
			CM_("ncrna", scolor::CYAN),
			CM_("variant", scolor::RED),
			CM_("mutation", scolor::RED),
			CM_("misc", scolor::VIOLET),
			CM_("misc+", scolor::DEEPSKY),
			CM_("misc-", scolor::CRIMSON)
		};
		//if (pref.hasKey("color-map")) {}		
		// Chromosome
		{
			SFigure chrArea;
			chrArea.setId("Chromosome");
			SLine chrBar(svec2f(0.f, CHR_BAR_OFFSET), svec2f(width, CHR_BAR_OFFSET), Stroke(sstyle::SOLID_LINE, 2.f));
			SCaption chrLabel(X_OFFSET, MAIN_TEXT_SIZE, adb.chromosomes[pos.idx].name);
			SFigure scaleLabels;
			Path2D<float> shortScale, longScale;
			//
			int lbin = smath::power(10, log(pos.length()) / log(10.0));
			if (pos.length() / lbin < 3) lbin /= 2;
			int unit = (1000000 <= lbin ? 1000000 : 1000), sbin = lbin / 10;
			String unitlabel = unit < 1000 ? "" : (unit < 10000 ? "K" : "M");
			//
			int lsc = (pos.begin / lbin + 1) * lbin, ssc = (pos.begin / sbin + 1) * sbin;
			//
			TextAttribute sclAttr;
			sclAttr.size = NOTE_TEXT_SIZE;
			//
			while (ssc < lsc) {
				shortScale.moveTo(svec2f((float)(ssc - pos.begin) * ratio, CHR_BAR_OFFSET - SSCALE_HEIGHT / 2.f));
				shortScale.lineTo(svec2f((float)(ssc - pos.begin) * ratio, CHR_BAR_OFFSET + SSCALE_HEIGHT / 2.f));
				ssc += sbin;
			}
			while (lsc < pos.end) {
				auto scLabel = SCaption((float)(lsc - pos.begin) * ratio - sclAttr.size / 2.f, CHR_BAR_OFFSET - LSCALE_HEIGHT / 2.f - NOTE_TEXT_SIZE / 2.f, S(lsc / unit) + unitlabel, sclAttr);
				scaleLabels.draw(scLabel);
				longScale.moveTo(svec2f((float)(lsc - pos.begin) * ratio, CHR_BAR_OFFSET - LSCALE_HEIGHT / 2.f));
				longScale.lineTo(svec2f((float)(lsc - pos.begin) * ratio, CHR_BAR_OFFSET + LSCALE_HEIGHT / 2.f));
				ssc = lsc + sbin;
				sforin(i, 1, 10) {
					if (pos.end < ssc) break;
					shortScale.moveTo(svec2f((float)(ssc - pos.begin) * ratio, CHR_BAR_OFFSET - SSCALE_HEIGHT / 2.f));
					shortScale.lineTo(svec2f((float)(ssc - pos.begin) * ratio, CHR_BAR_OFFSET + SSCALE_HEIGHT / 2.f));
					ssc += sbin;
				}
				lsc += lbin;
			}
			chrArea
				.draw(chrLabel)
				.draw(chrBar)
				.draw(scaleLabels)
				.draw(SPath(shortScale))
				.draw(SPath(longScale));
			gmap.draw(chrArea);
		}

		// Gene
		if (display & (sushort)ANNOT_CATEGORY::GENE || display & (sushort)ANNOT_CATEGORY::TRANSCRIPT) {
			auto& genes = adb.getGenes(pos, {D_("transcript", true)});
			// Display genes
			if (display & (sushort)ANNOT_CATEGORY::GENE) {
				SFigure geneArea;
				geneArea.setId("Gene");
				sforeach(gene, genes) {
					SFigure genefig,
						genelabel = SCaption((gene.begin < pos.begin ? X_OFFSET : (gene.begin - pos.begin) * ratio), Y_OFFSET, gene.name, boxAttr),
						genebox = SRectangle((gene.begin - pos.begin) * ratio, Y_OFFSET + LABEL_MARGIN, (float)gene.length(true) * ratio, BOX_HEIGHT);
					if (gene.type == (suint)GENE_TYPE::PSEUDO_GENE) genebox.setFillColor(colorMap["pseudo"]);
					else if (gene.type == (suint)GENE_TYPE::TRANSPOSON) genebox.setFillColor(colorMap["me"]);
					else genebox.setFillColor(colorMap["gene"]);
					genefig
						.draw(genelabel)
						.draw(genebox);
					genefig.setId(gene.geneid);
					geneArea.draw(genefig);
				}
				moir::rearrange(geneArea, VERTICAL_SHIFT);
				gmap.draw(geneArea);
				Y_OFFSET = gmap.boundary().ori_y + gmap.boundary().height + SECTOR_MARGIN;
			}
			// Display transcripts
			if (display & (sushort)ANNOT_CATEGORY::TRANSCRIPT){
				SFigure trsArea;
				trsArea.setId("Transcript");
				sforeach(gene, genes) {
					sforeach(trs, gene.transcripts) {
						SFigure trsfig,
							trslabel = SCaption((trs->begin < pos.begin ? 0 : (trs->begin - pos.begin)) * ratio + X_OFFSET, Y_OFFSET, trs->name, boxAttr),
							trsline = SLine((trs->begin - pos.begin) * ratio, Y_OFFSET + LABEL_MARGIN + BOX_HEIGHT / 2.f,
								(trs->end - pos.begin) * ratio, Y_OFFSET + LABEL_MARGIN + BOX_HEIGHT / 2.f),
							trsdir = SPolygon(3);
						SFigure trsboxs;
						sforeach(strct, trs->structures) {
							SRectangle strctbox((strct.begin - pos.begin) * ratio, Y_OFFSET + LABEL_MARGIN, (float)strct.length(true) * ratio, BOX_HEIGHT);
							if (gene.type == (suint)GENE_TYPE::PROTEIN_CODING) {
								if (strct.type == CDS) strctbox.brush = colorMap[(gene.dir ? "CDS-" : "CDS+")];
								else if (strct.type & UTR) strctbox.brush = colorMap["UTR"];
								else strctbox.brush = scolor::CLEAR;
							}
							else {
								if (strct.type == EXON) {
									if (trs->type == (suint)TRANSCRIPT_TYPE::T_RNA) strctbox.brush = colorMap["trna"];
									else if (trs->type == (suint)TRANSCRIPT_TYPE::R_RNA) strctbox.brush = colorMap["rrna"];
									else strctbox.brush = colorMap["ncrna"];
								}
							}
							trsboxs.draw(strctbox);
						}
						if (gene.dir) {
							trsdir.setVertex(0, svec2f((trs->begin - pos.begin)* ratio, Y_OFFSET + LABEL_MARGIN));
							trsdir.setVertex(1, svec2f((trs->begin - pos.begin)* ratio - arrowwidth, Y_OFFSET + LABEL_MARGIN + BOX_HEIGHT / 2));
							trsdir.setVertex(2, svec2f((trs->begin - pos.begin)* ratio, Y_OFFSET + LABEL_MARGIN + BOX_HEIGHT));
						}
						else {
							trsdir.setVertex(0, svec2f((trs->end - pos.begin)* ratio, Y_OFFSET + LABEL_MARGIN));
							trsdir.setVertex(1, svec2f((trs->end - pos.begin)* ratio, Y_OFFSET + LABEL_MARGIN + BOX_HEIGHT));
							trsdir.setVertex(2, svec2f((trs->end - pos.begin)* ratio + arrowwidth, Y_OFFSET + LABEL_MARGIN + BOX_HEIGHT / 2));
						}
						trsdir.setFillColor(scolor::BLACK);
						trsfig.setId(trs->name);
						trsfig
							.draw(trslabel)
							.draw(trsdir)
							.draw(trsline)
							.draw(trsboxs);
						trsArea.draw(trsfig);
					}
				}
				moir::rearrange(trsArea, VERTICAL_SHIFT);
				gmap.draw(trsArea);
				Y_OFFSET = gmap.boundary().ori_y + gmap.boundary().height + SECTOR_MARGIN;
			}
		}
		if (display & (sushort)ANNOT_CATEGORY::MUTATION) {
			SFigure varArea;
			varArea.setId("Mutation/Variant");
			auto vars = adb.getMutations(pos);
			vars.append(adb.getVariants(pos));
			sforeach(var, vars) {
				SFigure varfig,
					varlabel = SCaption((var.begin < pos.begin ? X_OFFSET : (var.begin - pos.begin) * ratio), Y_OFFSET, var.name, boxAttr),
					varbox = SRectangle((var.begin - pos.begin) * ratio, Y_OFFSET + LABEL_MARGIN, (float)var.length(true) * ratio, BOX_HEIGHT);
				varbox.setFillColor(colorMap["variant"]);
				varfig
					.draw(varlabel)
					.draw(varbox);
				varfig.setId(var.name);
				varArea.draw(varfig);
			}
			moir::rearrange(varArea, VERTICAL_SHIFT);
			gmap.draw(varArea);
			Y_OFFSET = gmap.boundary().ori_y + gmap.boundary().height + SECTOR_MARGIN;
		}
		if (display & (sushort)ANNOT_CATEGORY::FEATURE) {
			SFigure ftrArea;
			/*
			* 
			*/
		}
		if (display & (sushort)ANNOT_CATEGORY::USER_DEFINED) {
			SFigure customArea;
			auto objs = pref["mapping"].split(",");
			sforeach(args, objs) {
				auto obj = args.parse(";", "=");
				SFigure customFig;
				auto start = obj["pos"].intValue(), len = obj["length"].intValue();
				srange range = srange(start, start + len -1);
				int direction = obj.hasKey("dir") ? (obj["dir"] == "+" ? 0 : 1) : -1;

				SCaption customLabel((range.begin < pos.begin ? 0 : (range.begin - pos.begin)) * ratio, Y_OFFSET, obj["label"]);
				customFig.draw(customLabel);
				Color col = obj.hasKey("color") ? 
					Color(obj["color"]) : (direction == -1 ? colorMap["misc"]
						: (direction == 0 ? colorMap["misc+"] : colorMap["misc-"]));
				if (!obj.hasKey("type")) obj["type"] = "rect";
				if (obj["type"] == "line") {
					auto customLine = SLine(
						svec2f((range.begin - pos.begin) * ratio, Y_OFFSET + 15.f),
						svec2f((range.end - pos.begin) * ratio, Y_OFFSET + 15.f)
					);
					if (obj.hasKey("width")) customLine.stroke.width = obj["width"].floatValue();
					customLine.stroke.color = col;
					customFig.draw(customLine);
				}
				else if (obj["type"] == "rect") {
					SFigure customRect = SRectangle((range.begin - pos.begin) * ratio, Y_OFFSET + LABEL_MARGIN, range.length(true) * ratio, BOX_HEIGHT);
					customRect.setFillColor(col);
					if (-1 < direction) {
						SFigure customdir = SPolygon(3);
						if (direction == 0) {
							customdir.setVertex(0, svec2f((range.end - pos.begin) * ratio, Y_OFFSET + LABEL_MARGIN));
							customdir.setVertex(1, svec2f((range.end - pos.begin) * ratio, Y_OFFSET + LABEL_MARGIN + BOX_HEIGHT));
							customdir.setVertex(2, svec2f((range.end - pos.begin) * ratio + arrowwidth, Y_OFFSET + LABEL_MARGIN + BOX_HEIGHT / 2));
						}
						else {
							customdir.setVertex(0, svec2f((range.begin - pos.begin) * ratio, Y_OFFSET + LABEL_MARGIN));
							customdir.setVertex(1, svec2f((range.begin - pos.begin) * ratio - arrowwidth, Y_OFFSET + LABEL_MARGIN + BOX_HEIGHT / 2));
							customdir.setVertex(2, svec2f((range.begin - pos.begin) * ratio, Y_OFFSET + LABEL_MARGIN + BOX_HEIGHT));
						}
						customdir.setFillColor(scolor::BLACK);
					}
					customFig.draw(customRect);
				}
				customArea.draw(customFig);
			}
			moir::rearrange(customArea, VERTICAL_SHIFT);
			gmap.draw(customArea);
		}
		gmap.resize(width, gmap.boundary().height + SECTOR_MARGIN);
		if (param.oformat == "_obj_") param.response.attribute["result"] = gmap;
		else {
			if (param.output.empty()) {
				param.ostream.setStrOStream(param.response.output);
				param.ostream.print(sxml::svgNode(gmap)->toString());
			}
			else gmap.save(param.output);
			if (param.response.output && !pref["silent"]) SPrint(param.response.output);
		}
	}
	catch (Exception ex) {
		ex.print();
		param.response = ex;
	}
	return param.response;
}