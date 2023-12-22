/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2023 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

#include "internal/CustomSummationImpl.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/Platform.h"
#include "openmm/State.h"
#include "openmm/System.h"
#include "openmm/Vec3.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/ContextImpl.h"

#include <map>
#include <regex>
#include <string>

using namespace OpenMM;
using namespace OpenMMLab;
using namespace std;

CustomSummationImpl::CustomSummationImpl(
    int numArgs,
    CustomCompoundBondForce &force,
    Platform &platform,
    const map<string, string> &properties
) : numArgs(numArgs) {
    int numParticles = force.getNumParticlesPerBond();
    positions.resize(numParticles, Vec3(0, 0, 0));
    latestArguments.resize(numArgs);
    derivatives.resize(numArgs);

    System *system = new System();
    for (int i = 0; i < numParticles; i++)
        system->addParticle(1.0);
    system->addForce(static_cast<Force *>(&force));

    VerletIntegrator *integrator = new VerletIntegrator(0.01);
    context = new Context(*system, *integrator, platform, properties);
    valueIsDirty = derivativesAreDirty = contextIsUnchanged = true;
};

CustomSummationImpl::~CustomSummationImpl() {
    delete context; // will delete the system, integrator, and force
}

void CustomSummationImpl::setPositions(const vector<double> &arguments) {
    if (equal(arguments.begin(), arguments.end(), latestArguments.begin()) && contextIsUnchanged)
        return;
    for (int i = 0; i < numArgs; i++)
        positions[i / 3][i % 3] = arguments[i];
    context->setPositions(positions);
    latestArguments = arguments;
    valueIsDirty = derivativesAreDirty = contextIsUnchanged = true;
}

double CustomSummationImpl::evaluate(const vector<double> &arguments) {
    setPositions(arguments);
    if (valueIsDirty) {
        value = context->getState(State::Energy).getPotentialEnergy();
        valueIsDirty = false;
    }
    return value;
}

vector<double> CustomSummationImpl::evaluateDerivatives(const vector<double> &arguments) {
    setPositions(arguments);
    if (derivativesAreDirty) {
        vector<Vec3> forces = context->getState(State::Forces).getForces();
        for (int i = 0; i < numArgs; i++)
            derivatives[i] = -forces[i / 3][i % 3];
        derivativesAreDirty = false;
    }
    return derivatives;
}

void CustomSummationImpl::update(CustomCompoundBondForce &force) {
    force.updateParametersInContext(*context);
    contextIsUnchanged = false;
}

void CustomSummationImpl::reset() {
    context->reinitialize();
    contextIsUnchanged = false;
}

double CustomSummationImpl::getParameter(const string &name) const {
    return context->getParameter(name);
}

void CustomSummationImpl::setParameter(const string &name, double value) {
    context->setParameter(name, value);
    contextIsUnchanged = false;
}

const map<string, string> &CustomSummationImpl::getPlatformProperties() const {
    Platform &platform = context->getPlatform();
    vector<string> names = platform.getPropertyNames();
    map<string, string> &properties = *new map<string, string>;
    for (const auto& name : names)
        properties[name] = platform.getPropertyValue(*context, name);
    return properties;
}
