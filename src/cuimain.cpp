#include "sapp/scuiapp.h"
#include "profile.h"
#include "moirei.h"
class MoireiCUI : public slib::sapp::SCuiApp {
	Moirei moirei;
public:
	MoireiCUI() : slib::sapp::SCuiApp(app_profile, prof_format) {}
	~MoireiCUI() {}
	int exec() { return moirei.run(preference); }
};
RUN_CUI_APP(MoireiCUI)
