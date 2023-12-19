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
#include "openmm/Force.h"
#include "lepton/CustomFunction.h"
#include <string>
#include <vector>

using namespace OpenMM;
using namespace Lepton;
using namespace std;

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
        return arguments.size();
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
     * Add a new force to this function.
     *
     * @param force    the force to add
     *
     * @return the index of the force
     */
    int addForce(Force* force);
    /**
     * Get the force at the given index.
     *
     * @param index    the index of the force to get
     *
     * @return the force at the given index
     */
    Force& getForce(int index);
    /**
     * Add a new argument to this function.
     *
     * @param particle  the particle index associated with the argument
     * @param direction the spatial direction associated with the argument
     *
     * @return the index of the argument
     */
    int addArgument(int particle, const std::string& direction);
    /**
     * Get the argument at the given index.
     *
     * @param index    the index of the argument to get
     *
     * @return the argument at the given index
     */
    void getArgument(int index, int* particle, std::string* direction);
private:
    class ArgumentInfo;
    vector<Force*> forces;
    vector<ArgumentInfo> arguments;
};

} // namespace OpenMMLab

#endif /*OPENMMLAB_FORCEDERIVEDFUNCTION_H_*/
