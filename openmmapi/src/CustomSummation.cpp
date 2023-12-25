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
#include "openmm/OpenMMException.h"
#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"

#include <map>
#include <string>
#include <vector>

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
    impl = new CustomSummationImpl(
        numArgs, expression, overallParameters, perTermParameters, platform, platformProperties
    );
}

CustomSummation::~CustomSummation() {
    delete impl;
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
    CustomSummation *copy = new CustomSummation(
        numArgs, expression, overallParameters, perTermParameters, *platform, platformProperties
    );
    for (int i = 0; i < getNumTerms(); i++)
        copy->addTerm(getTerm(i));
    copy->reinitialize();
    return copy;
}

int CustomSummation::addTerm(const vector<double> &parameters) {
    ASSERT_EQUAL(parameters.size(), perTermParameters.size());
    termParameters.push_back(parameters);
    return termParameters.size() - 1;
}

const vector<double> &CustomSummation::getTerm(int index) const {
    ASSERT_VALID_INDEX(index, termParameters);
    return termParameters[index];
}

void CustomSummation::setTerm(int index, const vector<double> &parameters) {
    ASSERT_VALID_INDEX(index, termParameters);
    termParameters[index] = parameters;
}

double CustomSummation::getParameter(const string& name) const {
    auto it = overallParameters.find(name);
    if (it == overallParameters.end())
        throw OpenMMException("Unknown parameter '" + name + "'");
    return it->second;
}

void CustomSummation::setParameter(const string &name, double value) {
    auto it = overallParameters.find(name);
    if (it == overallParameters.end())
        throw OpenMMException("Unknown parameter '" + name + "'");
    it->second = value;
    impl->setParameter(name, value);
}

void CustomSummation::reinitialize() {
    impl->reset(termParameters);
}

void CustomSummation::update() {
    impl->reset(termParameters);
}
