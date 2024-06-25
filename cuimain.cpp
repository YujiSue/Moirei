#include "sapp/scuiapp.h"
#include "profile.h"
#include "moirei.h"
class MoireiCUI : public slib::sapp::SCuiApp {
	Moirei moirei;
public:
	MoireiCUI() : slib::sapp::SCuiApp(app_profile, prof_format) {}
	~MoireiCUI() {}
	int exec() {
		auto res = moirei.run(preference);
		if (!res.code) SPrint(res.output);
		return res.code;
	}
};
RUN_CUI_APP(MoireiCUI)
