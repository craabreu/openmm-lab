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
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/internal/AssertionUtilities.h"
#include <map>
#include <regex>
#include <string>

#define ASSERT_INDEX(index, num) { \
    if (index < 0 || index >= num) throwException(__FILE__, __LINE__, "Index out of range"); \
};

using namespace OpenMM;
using namespace OpenMMLab;
using namespace std;

string replaceSymbols(const string& expression, int numArgs, bool forward = true) {
    const string xyz = "xyz";
    string newExpression = expression;
    for (int i = 0; i < numArgs; ++i) {
        string oldSymbol, newSymbol;
        if (forward) {
            oldSymbol = "x" + to_string(i + 1);
            newSymbol = xyz[i % 3] + to_string(i / 3 + 1);
        }
        else {
            int j = numArgs - 1 - i;
            oldSymbol = xyz[j % 3] + to_string(j / 3 + 1);
            newSymbol = "x" + to_string(j + 1);
        }
        newExpression = regex_replace(
            newExpression, regex("\\b" + oldSymbol + "\\b"), newSymbol
        );
    }
    return newExpression;
}

CustomSummation::CustomSummation(
    int numArgs,
    const std::string &expression,
    const map<string, double> &overallParameters,
    const vector<string> &perTermParameters,
    Platform &platform,
    const map<string, string> &properties
) : numArgs(numArgs) {

    int numParticles = (numArgs + 2) / 3;
    particles = vector<int>(numParticles);
    System *system = new System();
    for (int i = 0; i < numParticles; i++) {
        particles[i] = i;
        system->addParticle(1.0);
    }

    force = new CustomCompoundBondForce(numParticles, replaceSymbols(expression, numArgs));
    force->setUsesPeriodicBoundaryConditions(false);
    for (const auto& pair : overallParameters)
        force->addGlobalParameter(pair.first, pair.second);
    for (const auto& name : perTermParameters)
        force->addPerBondParameter(name);
    system->addForce(force);

    VerletIntegrator *integrator = new VerletIntegrator(0.01);
    context = new Context(*system, *integrator, platform, properties);
}

CustomSummation::~CustomSummation() {
    delete force;
    delete context;
}

double CustomSummation::evaluate(const double* arguments) const {
    return 0.0;
}

double CustomSummation::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    return 0.0;
}

CustomSummation* CustomSummation::clone() const {
    Platform *platform = &context->getPlatform();
    string expression = replaceSymbols(force->getEnergyFunction(), numArgs, false);
    map<string, double> overallParameters;
    for (int i = 0; i < force->getNumGlobalParameters(); i++)
        overallParameters[
            force->getGlobalParameterName(i)
        ] = force->getGlobalParameterDefaultValue(i);
    vector<string> perTermParameters;
    for (int i = 0; i < force->getNumPerBondParameters(); i++)
        perTermParameters.push_back(force->getPerBondParameterName(i));
    CustomSummation *summation = new CustomSummation(
        numArgs, expression, overallParameters, perTermParameters, *platform
    );
    for (int i = 0; i < force->getNumGlobalParameters(); i++)
        summation->context->setParameter(
            force->getGlobalParameterName(i),
            context->getParameter(force->getGlobalParameterName(i))
        );
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
    return context->getParameter(force->getGlobalParameterName(index));
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
    context->reinitialize();
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
    force->updateParametersInContext(*context);
}
