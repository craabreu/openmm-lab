/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2023 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

#include "ForceDerivedFunction.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"

using namespace OpenMMLab;

class ForceDerivedFunction::ArgumentInfo {
public:
    int particle;
    std::string direction;
    ArgumentInfo() {}
    ArgumentInfo(int particle, const std::string& direction) : particle(particle), direction(direction) {}
};

ForceDerivedFunction::ForceDerivedFunction() {
    forces.resize(0);
    arguments.resize(0);
}

ForceDerivedFunction::~ForceDerivedFunction() {
}

double ForceDerivedFunction::evaluate(const double* arguments) const {
    return 0.0;
}

double ForceDerivedFunction::evaluateDerivative(const double* arguments, const int* derivOrder) {
    return 0.0;
}

int ForceDerivedFunction::addForce(Force* force) {
    forces.push_back(force);
    return forces.size()-1;
}

Force& ForceDerivedFunction::getForce(int index) {
    ASSERT_VALID_INDEX(index, forces);
    return *forces[index];
}

int ForceDerivedFunction::addArgument(int particle, const std::string& direction) {
    arguments.push_back(ArgumentInfo(particle, direction));
    return arguments.size()-1;
}

void ForceDerivedFunction::getArgument(int index, int* particle, std::string* direction) {
    ASSERT_VALID_INDEX(index, arguments);
    *particle = arguments[index].particle;
    *direction = arguments[index].direction;
}
