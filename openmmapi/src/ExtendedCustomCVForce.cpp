/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2023 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

#include "ExtendedCustomCVForce.h"
#include "internal/ExtendedCustomCVForceImpl.h"
#include "OpenMMLabKernels.h"

#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include <cmath>
#include <map>
#include <set>
#include <utility>

using namespace OpenMMLab;
using namespace OpenMM;
using namespace std;

ExtendedCustomCVForce::ExtendedCustomCVForce(const string& energy) : energyExpression(energy) {
    this->setName("ExtendedCustomCVForce");
}

ExtendedCustomCVForce::~ExtendedCustomCVForce() {
    for (auto variable : variables)
        delete variable.variable;
    for (auto function : functions)
        delete function.function;
}

const string& ExtendedCustomCVForce::getEnergyFunction() const {
    return energyExpression;
}

void ExtendedCustomCVForce::setEnergyFunction(const std::string& energy) {
    energyExpression = energy;
}

int ExtendedCustomCVForce::addCollectiveVariable(const std::string& name, Force* variable) {
    if (variables.size() >= 32)
        throw OpenMMException("ExtendedCustomCVForce cannot have more than 32 collective variables");
    variables.push_back(VariableInfo(name, variable));
    return variables.size()-1;
}

const string& ExtendedCustomCVForce::getCollectiveVariableName(int index) const {
    ASSERT_VALID_INDEX(index, variables);
    return variables[index].name;
}

Force& ExtendedCustomCVForce::getCollectiveVariable(int index) {
    ASSERT_VALID_INDEX(index, variables);
    return *variables[index].variable;
}

const Force& ExtendedCustomCVForce::getCollectiveVariable(int index) const {
    ASSERT_VALID_INDEX(index, variables);
    return *variables[index].variable;
}

int ExtendedCustomCVForce::addGlobalParameter(const string& name, double defaultValue) {
    globalParameters.push_back(GlobalParameterInfo(name, defaultValue));
    return globalParameters.size()-1;
}

const string& ExtendedCustomCVForce::getGlobalParameterName(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].name;
}

void ExtendedCustomCVForce::setGlobalParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].name = name;
}

double ExtendedCustomCVForce::getGlobalParameterDefaultValue(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].defaultValue;
}

void ExtendedCustomCVForce::setGlobalParameterDefaultValue(int index, double defaultValue) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].defaultValue = defaultValue;
}

void ExtendedCustomCVForce::addEnergyParameterDerivative(const string& name) {
    for (int i = 0; i < globalParameters.size(); i++)
        if (name == globalParameters[i].name) {
            energyParameterDerivatives.push_back(i);
            return;
        }
    throw OpenMMException(string("addEnergyParameterDerivative: Unknown global parameter '"+name+"'"));
}

const string& ExtendedCustomCVForce::getEnergyParameterDerivativeName(int index) const {
    ASSERT_VALID_INDEX(index, energyParameterDerivatives);
    return globalParameters[energyParameterDerivatives[index]].name;
}

int ExtendedCustomCVForce::addTabulatedFunction(const std::string& name, TabulatedFunction* function) {
    functions.push_back(FunctionInfo(name, function));
    return functions.size()-1;
}

const TabulatedFunction& ExtendedCustomCVForce::getTabulatedFunction(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

TabulatedFunction& ExtendedCustomCVForce::getTabulatedFunction(int index) {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

const string& ExtendedCustomCVForce::getTabulatedFunctionName(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return functions[index].name;
}

void ExtendedCustomCVForce::getCollectiveVariableValues(Context& context, vector<double>& values) const {
    dynamic_cast<ExtendedCustomCVForceImpl&>(getImplInContext(context)).getCollectiveVariableValues(getContextImpl(context), values);
}

ForceImpl* ExtendedCustomCVForce::createImpl() const {
    return new ExtendedCustomCVForceImpl(*this);
}

Context& ExtendedCustomCVForce::getInnerContext(Context& context) {
    return dynamic_cast<ExtendedCustomCVForceImpl&>(getImplInContext(context)).getInnerContext();
}

void ExtendedCustomCVForce::updateParametersInContext(Context& context) {
    dynamic_cast<ExtendedCustomCVForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

bool ExtendedCustomCVForce::usesPeriodicBoundaryConditions() const {
    for (auto& variable : variables)
        if (variable.variable->usesPeriodicBoundaryConditions())
            return true;
    return false;
}
