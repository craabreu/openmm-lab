/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2023 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

#include "CustomSummation.h"
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

#include <iostream>  // temporary

#define ASSERT_INDEX(index, num) { \
    if (index < 0 || index >= num) throwException(__FILE__, __LINE__, "Index out of range"); \
};

using namespace OpenMM;
using namespace OpenMMLab;
using namespace std;

CustomSummation::Evaluator::Evaluator(
    int numArgs,
    CustomCompoundBondForce &force,
    Platform &platform,
    const map<string, string> &properties
) : numArgs(numArgs), valueIsDirty(true), derivativesAreDirty(true) {
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
};

CustomSummation::Evaluator::~Evaluator() {
    delete context; // will delete the system, integrator, and force
}

void CustomSummation::Evaluator::setPositions(vector<double> arguments) {
    if (equal(arguments.begin(), arguments.end(), latestArguments.begin()))
        return;
    for (int i = 0; i < numArgs; i++)
        positions[i / 3][i % 3] = arguments[i];
    context->setPositions(positions);
    latestArguments = arguments;
    valueIsDirty = true;
    derivativesAreDirty = true;
}

double CustomSummation::Evaluator::evaluate(vector<double> arguments) {
    setPositions(arguments);
    if (valueIsDirty) {
        value = context->getState(State::Energy).getPotentialEnergy();
        valueIsDirty = false;
    }
    return value;
}

vector<double> CustomSummation::Evaluator::evaluateDerivatives(vector<double> arguments) {
    setPositions(arguments);
    if (derivativesAreDirty) {
        vector<Vec3> forces = context->getState(State::Forces).getForces();
        for (int i = 0; i < numArgs; i++)
            derivatives[i] = -forces[i / 3][i % 3];
        derivativesAreDirty = false;
    }
    return derivatives;
}

CustomSummation::CustomSummation(
    int numArgs,
    const std::string &expression,
    const map<string, double> &overallParameters,
    const vector<string> &perTermParameters,
    Platform &platform,
    const map<string, string> &properties
) : numArgs(numArgs) {
    int numParticles = (numArgs  + 2)/ 3;
    for (int i = 0; i < numParticles; i++)
        particles.push_back(i);
    force = new CustomCompoundBondForce(numParticles, expression);
    force->setUsesPeriodicBoundaryConditions(false);
    for (const auto& pair : overallParameters)
        force->addGlobalParameter(pair.first, pair.second);
    for (const auto& name : perTermParameters)
        force->addPerBondParameter(name);
    evaluator = new Evaluator(numArgs, *force, platform, properties);
}

CustomSummation::~CustomSummation() {
    delete evaluator;  // will delete the force
}

double CustomSummation::evaluate(const double* arguments) const {
    return evaluator->evaluate(vector<double>(arguments, arguments + numArgs));
}

double CustomSummation::evaluate(const vector<double> &arguments) const {
    return evaluator->evaluate(arguments);
}

double CustomSummation::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    int index;
    int sum = 0;
    for (int i = 0; i < numArgs; i++) {
        if (derivOrder[i] < 0 || derivOrder[i] > 1)
            throwException(__FILE__, __LINE__, "Invalid derivative order");
        if (derivOrder[i] == 1) {
            index = i;
            sum++;
        }
    }
    if (sum != 1)
        throwException(__FILE__, __LINE__, "Exactly one derivative must be specified");
    vector<double> derivatives = evaluator->evaluateDerivatives(
        vector<double>(arguments, arguments + numArgs)
    );
    return derivatives[index];
}

double CustomSummation::evaluateDerivative(const vector<double> &arguments, int which) const {
    return evaluator->evaluateDerivatives(arguments)[which];
}

CustomSummation* CustomSummation::clone() const {
    map<string, double> overallParameters;
    for (int i = 0; i < force->getNumGlobalParameters(); i++)
        overallParameters[
            force->getGlobalParameterName(i)
        ] = force->getGlobalParameterDefaultValue(i);
    vector<string> perTermParameters;
    for (int i = 0; i < force->getNumPerBondParameters(); i++)
        perTermParameters.push_back(force->getPerBondParameterName(i));
    CustomSummation *summation = new CustomSummation(
        numArgs,
        force->getEnergyFunction(),
        overallParameters,
        perTermParameters,
        evaluator->context->getPlatform()
    );
    for (int i = 0; i < force->getNumGlobalParameters(); i++) {
        string name = force->getGlobalParameterName(i);
        summation->setOverallParameter(name, evaluator->context->getParameter(name));
    }
    return summation;
}

int CustomSummation::getNumOverallParameters() const {
    return force->getNumGlobalParameters();
}

const string& CustomSummation::getOverallParameterName(int index) const {
    ASSERT_INDEX(index, force->getNumGlobalParameters());
    return force->getGlobalParameterName(index);
}

double CustomSummation::getOverallParameterDefaultValue(int index) const {
    ASSERT_INDEX(index, force->getNumGlobalParameters());
    return evaluator->context->getParameter(force->getGlobalParameterName(index));
}

int CustomSummation::getNumPerTermParameters() const {
    return force->getNumPerBondParameters();
}

const string& CustomSummation::getPerTermParameterName(int index) const {
    ASSERT_INDEX(index, force->getNumPerBondParameters());
    return force->getPerBondParameterName(index);
}

int CustomSummation::addTerm(vector<double> parameters) {
    force->addBond(particles, parameters);
    evaluator->context->reinitialize();
    return force->getNumBonds() - 1;
}

vector<double> CustomSummation::getTermParameters(int index) const {
    ASSERT_INDEX(index, force->getNumBonds());
    vector<int> particles;
    vector<double> parameters;
    force->getBondParameters(index, particles, parameters);
    return parameters;
}

void CustomSummation::setTermParameters(int index, vector<double> parameters) {
    ASSERT_INDEX(index, force->getNumBonds());
    force->setBondParameters(index, particles, parameters);
    force->updateParametersInContext(*evaluator->context);
}

void CustomSummation::setOverallParameter(const string &name, double value) {
    evaluator->context->setParameter(name, value);
}
