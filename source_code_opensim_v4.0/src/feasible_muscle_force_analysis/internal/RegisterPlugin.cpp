#include "RegisterPlugin.h"

#include <OpenSim/Common/Object.h>
#include "../FeasibleMuscleForceAnalysis.h"

using namespace OpenSim;

static dllObjectInstantiator instantiator;

void RegisterPlugin() {
    Object::RegisterType(FeasibleMuscleForceAnalysis());
}

dllObjectInstantiator::dllObjectInstantiator() {
    registerDllClasses();
}

void dllObjectInstantiator::registerDllClasses() {
    RegisterPlugin();
}
