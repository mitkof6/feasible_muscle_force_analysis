#ifndef REGISTER_PLUGIN_H
#define REGISTER_PLUGIN_H

#include "FeasibleMuscleForceAnalysisExports.h"

extern "C" {
    /**
    * The purpose of this routine is to register all class types exported by
    * the plugin library.
    */
    FeasibleMuscleForceAnalysis_API void RegisterPlugin();
}

class dllObjectInstantiator {
public:
    dllObjectInstantiator();
private:
    void registerDllClasses();
};

#endif
