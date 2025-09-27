#ifndef _XRAY_HOST_CONTEXT
#define _XRAY_HOST_CONTEXT

#include "base_xrayHostContext.h"

typedef base_xrayHostContext* xrayHostContext;
typedef base_xrayHostContext _xrayHostContext;

// Singleton pattern for the global instance of an X-ray simulation class
class theXrayHostContext {

public:
  static _xrayHostContext& GetReference() {
    static _xrayHostContext thisInstance;
    static bool init = false;
    if (!init)  {
      thisInstance.Init();
      init = true;
    }
    return thisInstance;
  };

  static xrayHostContext GetPointer() {
    return (&GetReference());
  };

protected:
  theXrayHostContext(void);
  virtual ~theXrayHostContext(void);
};

#endif
