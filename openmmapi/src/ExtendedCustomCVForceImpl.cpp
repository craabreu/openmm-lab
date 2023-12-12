/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2023 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */


#include "internal/ExtendedCustomCVForceImpl.h"
#include "OpenMMLabKernels.h"

#include "openmm/NonbondedForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/serialization/XmlSerializer.h"
#include <map>

using namespace OpenMMLab;
using namespace OpenMM;
using namespace std;

ExtendedCustomCVForceImpl::ExtendedCustomCVForceImpl(const ExtendedCustomCVForce& owner) : owner(owner), innerIntegrator(1.0),
        innerContext(NULL) {
#if (OPENMM_VERSION_MAJOR > 8 || (OPENMM_VERSION_MAJOR == 8 && OPENMM_VERSION_MINOR > 0))
    forceGroup = owner.getForceGroup();
#endif
}

ExtendedCustomCVForceImpl::~ExtendedCustomCVForceImpl() {
    if (innerContext != NULL)
        delete innerContext;
}

void ExtendedCustomCVForceImpl::initialize(ContextImpl& context) {
    // Construct the inner system used to evaluate collective variables.

    const System& system = context.getSystem();
    Vec3 a, b, c;
    system.getDefaultPeriodicBoxVectors(a, b, c);
    innerSystem.setDefaultPeriodicBoxVectors(a, b, c);
    for (int i = 0; i < system.getNumParticles(); i++)
        innerSystem.addParticle(system.getParticleMass(i));
    for (int i = 0; i < owner.getNumCollectiveVariables(); i++) {
        Force* variable = XmlSerializer::clone<Force>(owner.getCollectiveVariable(i));
        variable->setForceGroup(i);
        NonbondedForce* nonbonded = dynamic_cast<NonbondedForce*>(variable);
        if (nonbonded != NULL)
            nonbonded->setReciprocalSpaceForceGroup(-1);
        innerSystem.addForce(variable);
    }

    // Create the inner context.

    innerContext = context.createLinkedContext(innerSystem, innerIntegrator);
    vector<Vec3> positions(system.getNumParticles(), Vec3());
    innerContext->setPositions(positions);

    // Create the kernel.

    kernel = context.getPlatform().createKernel(CalcExtendedCustomCVForceKernel::Name(), context);
    kernel.getAs<CalcExtendedCustomCVForceKernel>().initialize(context.getSystem(), owner, getContextImpl(*innerContext));
}

double ExtendedCustomCVForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
#if (OPENMM_VERSION_MAJOR > 8 || (OPENMM_VERSION_MAJOR == 8 && OPENMM_VERSION_MINOR > 0))
    if ((groups&(1<<forceGroup)) != 0)
#else
    if ((groups&(1<<owner.getForceGroup())) != 0)
#endif
        return kernel.getAs<CalcExtendedCustomCVForceKernel>().execute(context, getContextImpl(*innerContext), includeForces, includeEnergy);
    return 0.0;
}

vector<string> ExtendedCustomCVForceImpl::getKernelNames() {
    vector<string> names;
    names.push_back(CalcExtendedCustomCVForceKernel::Name());
    return names;
}

#if OPENMM_VERSION_MAJOR >= 8
vector<pair<int, int> > ExtendedCustomCVForceImpl::getBondedParticles() const {
    vector<pair<int, int> > bonds;
    const ContextImpl& innerContextImpl = getContextImpl(*innerContext);
    for (auto& impl : innerContextImpl.getForceImpls()) {
        for (auto& bond : impl->getBondedParticles())
            bonds.push_back(bond);
    }
    return bonds;
}
#endif

map<string, double> ExtendedCustomCVForceImpl::getDefaultParameters() {
    map<string, double> parameters;
    parameters.insert(innerContext->getParameters().begin(), innerContext->getParameters().end());
    for (int i = 0; i < owner.getNumGlobalParameters(); i++)
        parameters[owner.getGlobalParameterName(i)] = owner.getGlobalParameterDefaultValue(i);
    return parameters;
}

void ExtendedCustomCVForceImpl::getCollectiveVariableValues(ContextImpl& context, vector<double>& values) {
    kernel.getAs<CalcExtendedCustomCVForceKernel>().copyState(context, getContextImpl(*innerContext));
    values.clear();
    for (int i = 0; i < innerSystem.getNumForces(); i++) {
        double value = innerContext->getState(State::Energy, false, 1<<i).getPotentialEnergy();
        values.push_back(value);
    }
}

Context& ExtendedCustomCVForceImpl::getInnerContext() {
    return *innerContext;
}

void ExtendedCustomCVForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcExtendedCustomCVForceKernel>().copyParametersToContext(context, owner);
    context.systemChanged();
}
