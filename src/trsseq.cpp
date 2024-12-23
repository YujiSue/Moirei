#include "moirei.h"
inline void getTranscriptSeq(moir::Param& param, SeqList& reference, TranscriptInfo& transcript, SeqList& result, bool cds, bool rna) {
	Sequence seq;
	if (cds) {
		if (transcript.type == (int)TRANSCRIPT_TYPE::M_RNA)
			seq.setSeqAs(reference[transcript.gene->idx].raw(shift(transcript.coding(), -1), transcript.gene->dir), DNA_SEQ);
		else return;
	}
	else seq.setSeqAs(reference[transcript.gene->idx].raw(shift(transcript.exons(), -1), transcript.gene->dir), DNA_SEQ);
	if (rna) seq.transcribe();
	seq.name = transcript.name;
	//
	result.add(seq);
}
Response& Moirei::transcriptSeq(const SDictionary& pref) {
	try {
		param.setPref(pref);
		//
		SeqList reference, seqs;
		if (pref["verbose"]) SPrint("Load reference.");
		reference.load(pref["reference"]);
		//
		if (pref["verbose"]) SPrint("Connect annotation db.");
		// Annotation
		sushort atypes = 0;
		AnnotDB db;
		if (pref.hasKey("annotation") && pref["annotation"].size()) {
			if (pref["verbose"]) SPrint("Connect annotation database.");
			if (!pref.hasKey("annotdb")) throw InsufficientArgsException(insufficientArgsErrTxt("annotdb"));
			//
			auto str = sstr::toUpper(pref["annotation"]);
			sforeach(c, str) {
				switch (c) {
				case 'T':
					atypes |= (sushort)ANNOT_CATEGORY::TRANSCRIPT;
					break;
				case 'V':
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
		//
		db.open(pref["annotdb"]);
		bool rna = pref.hasKey("base-type") ? (pref["base-type"] == "rna") : false;
		bool gene = pref.hasKey("query-type")? (pref["query-type"] == "gene") : false;
		bool cds = pref.hasKey("cds-only") ? (bool)pref["cds-only"] : false;
		//
		sfor(pref["_args_"]) {
			if (gene) {
				auto& genes = db.getGenes($_, MATCH::PARTIAL);
				sforeach(gene, genes) {
					sforeach(transcript, gene.transcripts) {
						getTranscriptSeq(param, reference, *transcript, seqs, cds, rna);
					}
				}
			}
			else {
				auto& transcripts = db.getTranscripts($_, MATCH::PARTIAL);
				sforeach(transcript, transcripts) {
					getTranscriptSeq(param, reference, transcript, seqs, cds, rna);
				}
			}
		}
		if (pref["verbose"]) SPrint("Completed.");
		// Output
		if (param.oformat == "auto") param.oformat = "fa";
		moir::exportSeqs(param, seqs);
		if (param.response.output) SPrint(param.response.output);
	}
	catch (Exception ex) {
		ex.print();
		param.response = ex;
	}
	return param.response;
}
