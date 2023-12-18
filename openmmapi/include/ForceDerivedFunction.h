#ifndef OPENMMLAB_FORCEDERIVEDFUNCTION_H_
#define OPENMMLAB_FORCEDERIVEDFUNCTION_H_

/* -------------------------------------------------------------------------- *
 *                             OpenMM Laboratory                              *
 *                             =================                              *
 *                                                                            *
 * A plugin for testing low-level code implementation for OpenMM.             *
 *                                                                            *
 * Copyright (c) 2023 Charlles Abreu                                          *
 * https://github.com/craabreu/openmm-lab                                     *
 * -------------------------------------------------------------------------- */

#include "internal/windowsExportOpenMMLab.h"
#include "lepton/CustomFunction.h"

using namespace Lepton;

namespace OpenMMLab {

/*
 * A ForceDerivedFunction...
 */

class OPENMM_EXPORT_OPENMM_LAB ForceDerivedFunction : public Lepton::CustomFunction {
public:
    ForceDerivedFunction();
    ~ForceDerivedFunction();
    /**
     * Get the number of arguments this function expects.
     */
    int getNumArguments() const {
        return numArgs;
    };
    /**
     * Evaluate the function.
     *
     * @param arguments    the array of argument values
     */
    double evaluate(const double* arguments) const;
    /**
     * Evaluate a derivative of the function.
     *
     * @param arguments    the array of argument values
     * @param derivOrder   an array specifying the number of times the function has been differentiated
     *                     with respect to each of its arguments.  For example, the array {0, 2} indicates
     *                     a second derivative with respect to the second argument.
     */
    double evaluateDerivative(const double* arguments, const int* derivOrder);
    /**
     * Create a new duplicate of this object on the heap using the "new" operator.
     */
private:
    int numArgs;
};

} // namespace OpenMMLab

#endif /*OPENMMLAB_FORCEDERIVEDFUNCTION_H_*/
