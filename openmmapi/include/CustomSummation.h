#ifndef OPENMMLAB_CUSTOMSUMMATION_H_
#define OPENMMLAB_CUSTOMSUMMATION_H_

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
#include "openmm/Context.h"
#include "openmm/CustomExternalForce.h"
#include "openmm/Platform.h"
#include "lepton/CustomFunction.h"
#include <map>
#include <string>
#include <vector>

using namespace OpenMM;
using namespace Lepton;
using namespace std;

namespace OpenMMLab {

/*
* A CustomSummation...
*/

class OPENMM_EXPORT_OPENMM_LAB CustomSummation : public Lepton::CustomFunction {
public:
    /**
     * Construct a new CustomSummation object.
     *
     * @param numArgs             the number of arguments
     * @param expression          the expression for the function
     * @param overallParameters   the names and default values of the parameters that
     *                            are shared by all terms of the summation
     * @param perTermParameters   the names of the parameters that are unique to each
     *                            term of the summation
     * @param platform            the platform that will be used to evaluate the
     *                            summation
     * @param properties          a set of values for platform-specific properties.
     */
    CustomSummation(
        int numArgs,
        const std::string &expression,
        const map<string, double> &overallParameters,
        const vector<string> &perTermParameters,
        Platform &platform,
        const map<string, string> &properties = map<string, string>()
    );
    ~CustomSummation();
    /**
     * Get the number of arguments this function expects.
     */
    int getNumArguments() const { return numArgs; }
    /**
     * Evaluate the function.
     *
     * @param arguments    the array of argument values
     */
    double evaluate(const double *arguments) const;
    /**
     * Evaluate a derivative of the function.
     *
     * @param arguments    the array of argument values
     * @param derivOrder   an array specifying the number of times the function has been differentiated
     *                     with respect to each of its arguments.  For example, the array {0, 2} indicates
     *                     a second derivative with respect to the second argument.
     */
    double evaluateDerivative(const double *arguments, const int *derivOrder) const;
    /**
     * Create a new duplicate of this object on the heap using the "new" operator.
     */
    CustomSummation *clone() const;
    /**
     * Get the name of a overall parameter.
     *
     * @param index    the index of the overall parameter for which to get the name
     * @return         the overall parameter name
     */
    const std::string &getOverallParameterName(int index) const;
    /**
     * Get the value of a overall parameter.
     *
     * @param index    the index of the overall parameter for which to get the value
     * @return         the overall parameter value
     */
    double getOverallParameterDefaultValue(int index) const;
private:
    int numArgs;
    CustomExternalForce *force;
    Context *context;
};

} // namespace OpenMMLab

#endif /*OPENMMLAB_CUSTOMSUMMATION_H_*/
