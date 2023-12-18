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

using namespace OpenMMLab;

ForceDerivedFunction::ForceDerivedFunction() : numArgs(0) {
}

ForceDerivedFunction::~ForceDerivedFunction() {
}

double ForceDerivedFunction::evaluate(const double* arguments) const {
    return 0.0;
}

double ForceDerivedFunction::evaluateDerivative(const double* arguments, const int* derivOrder) {
    return 0.0;
}
