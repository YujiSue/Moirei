#include "moirei.h"
Response& Moirei::getComplement(const SDictionary& pref) {
	try {
		param.setPref(pref);
		//
		Sequence seq;
		seq.setSeq(pref["_args_"][0]);
		if (seq.type() != sseq::DNA) throw FormatException(formatErrorText("Sequence type", sbio::sutil::seqTypeName(seq.type()), "DNA"));
		//
		srange range(1, 0x7FFFFFFF);
		if (pref.hasKey("limit-query-length")) range.end = pref["limit-query-length"].intValue();
		if (seq.length() < range.begin || range.end < seq.length()) throw RangeException(outRangeErrorText("Sequence length", seq.length(), range.begin, range.end));
		//
		seq.complement();
		param.response.output = seq.raw();
		if (!pref["silent"]) SPrint(param.response.output);
	}
	catch (Exception ex) {
		ex.print();
		param.response = ex;
	}
	return param.response;
}
inline void loadCodon(CODON_TABLE& codon, const char* val) {
	String codons(val);
	if (codons == "default") codon = DEFAULT_CODON;
	else if (codons == "mt") codon = DEFAULT_MT_CODON;
	else {
		/*
		*/
	}
}
Response& Moirei::getTranslated(const SDictionary& pref) {
	try {
		param.setPref(pref);
		//
		CODON_TABLE codon;
		loadCodon(codon, pref["codon"]);
		//
		Sequence seq;
		seq.setSeq(pref["_args_"][0]);
		if (seq.type() != sseq::DNA) throw FormatException(formatErrorText("Sequence type", sbio::sutil::seqTypeName(seq.type()), "DNA"));
		//
		srange range(1, 0x7FFFFFFF);
		if (pref.hasKey("limit-query-length")) range.end = pref["limit-query-length"].intValue();
		if (seq.length()  < range.begin|| range.end < seq.length()) throw RangeException(outRangeErrorText("Sequence length", seq.length(), range.begin, range.end));
		//
		seq.transcribe().translate(codon);
		param.response.output = seq.raw();
		if (!pref["silent"]) SPrint(param.response.output);
	}
	catch (Exception ex) {
		ex.print();
		param.response = ex;
	}
	return param.response;
}
