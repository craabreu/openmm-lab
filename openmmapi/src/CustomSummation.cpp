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
    if (index < 0 || index >= num) \
        throwException(__FILE__, __LINE__, "Index out of range"); \
};

using namespace OpenMM;
using namespace OpenMMLab;
using namespace std;

CustomSummation::CustomSummation(
    int numArgs,
    string expression,
    map<string, double> overallParameters,
    vector<string> perTermParameters,
    Platform &platform,
    map<string, string> platformProperties
) : numArgs(numArgs),
    expression(expression),
    overallParameters(overallParameters),
    perTermParameters(perTermParameters),
    platform(&platform),
    platformProperties(platformProperties)
{
    int numParticles = (numArgs  + 2)/ 3;
    for (int i = 0; i < numParticles; i++)
        particles.push_back(i);
    force = new CustomCompoundBondForce(numParticles, expression);
    force->setUsesPeriodicBoundaryConditions(false);
    for (const auto& pair : overallParameters)
        force->addGlobalParameter(pair.first, pair.second);
    for (const auto& name : perTermParameters)
        force->addPerBondParameter(name);
    impl = new CustomSummationImpl(numArgs, *force, platform, platformProperties);
}

CustomSummation::~CustomSummation() {
    delete impl;  // will delete the force
}

double CustomSummation::evaluate(const double* arguments) const {
    return impl->evaluate(vector<double>(arguments, arguments + numArgs));
}

double CustomSummation::evaluate(const vector<double> &arguments) const {
    return impl->evaluate(arguments);
}

double CustomSummation::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    vector<int> order(derivOrder, derivOrder + numArgs);
    int which;
    int acc = 0;
    for (int i = 0; i < numArgs; i++) {
        int item = order[i];
        if (item < 0 || (acc += item) > 1)
            throwException(__FILE__, __LINE__, "Invalid derivative order specification");
        if (item == 1)
            which = i;
    }
    vector<double> args(arguments, arguments + numArgs);
    return impl->evaluateDerivatives(args)[which];
}

double CustomSummation::evaluateDerivative(const vector<double> &arguments, int which) const {
    return impl->evaluateDerivatives(arguments)[which];
}

CustomSummation* CustomSummation::clone() const {
    map<string, double> overallParameters;
    for (int i = 0; i < getNumOverallParameters(); i++) {
        string name = getOverallParameterName(i);
        overallParameters[name] = getOverallParameterDefaultValue(i);
    }
    vector<string> perTermParameters;
    for (int i = 0; i < getNumPerTermParameters(); i++)
        perTermParameters.push_back(getPerTermParameterName(i));
    CustomSummation *copy = new CustomSummation(
        numArgs,
        expression,
        overallParameters,
        perTermParameters,
        *platform,
        platformProperties
    );
    for (int i = 0; i < getNumTerms(); i++)
        copy->addTerm(getTerm(i));
    for (int i = 0; i < getNumOverallParameters(); i++) {
        string name = getOverallParameterName(i);
        copy->setParameter(name, getParameter(name));
    }
    return copy;
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
    return impl->getParameter(force->getGlobalParameterName(index));
}

int CustomSummation::getNumPerTermParameters() const {
    return force->getNumPerBondParameters();
}

const string& CustomSummation::getPerTermParameterName(int index) const {
    ASSERT_INDEX(index, force->getNumPerBondParameters());
    return force->getPerBondParameterName(index);
}

int CustomSummation::addTerm(const vector<double> &parameters) {
    force->addBond(particles, parameters);
    impl->reset();
    return force->getNumBonds() - 1;
}

const vector<double> &CustomSummation::getTerm(int index) const {
    ASSERT_INDEX(index, force->getNumBonds());
    vector<int> _;
    vector<double> *parameters = new vector<double>(getNumPerTermParameters());
    force->getBondParameters(index, _, *parameters);
    return *parameters;
}

void CustomSummation::setTerm(int index, const vector<double> &parameters) {
    ASSERT_INDEX(index, force->getNumBonds());
    force->setBondParameters(index, particles, parameters);
    impl->update(*force);
}

double CustomSummation::getParameter(const string &name) const {
    return impl->getParameter(name);
}

void CustomSummation::setParameter(const string &name, double value) {
    impl->setParameter(name, value);
}
