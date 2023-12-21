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
#include "openmm/CustomCompoundBondForce.h"
#include "openmm/Platform.h"
#include "lepton/CustomFunction.h"
#include "openmm/internal/ContextImpl.h"
#include <map>
#include <string>
#include <vector>

using namespace OpenMM;
using namespace Lepton;
using namespace std;

namespace OpenMMLab {

/**
 * This class allows users to define a custom function that can be evaluated by means of
 * an OpenMM::Platform. It defines a sum that depends of a fixed number of arguments, a
 * set of per-term parameters, and a set of overall parameters.
 *
 * We refer to the arguments as x1, y1, z1, x2, y2, z2, x3, y3, etc. CustomSummation
 * evaluates a user supplied algebraic expression to determine the value of each term.
 * The expression may depend on the following variables and functions:
 *
 * <ul>
 * <li>x1, y1, z1, x2, y2, z2, x3, etc.: the argument passed to the function.</li>
 * <li>p1, p2, p3, etc.: three-dimensional points defined as (x1, y1, z1), (x2, y2, z2),
 * etc. If the number of arguments is not a multiple of 3, the last point is completed
 * with zeros.</li>
 * <li>distance(p1, p2): the distance between points p1 and p2 (where "p1" and "p2"
 * may be replaced by any valid point names.</li>
 * <li>angle(p1, p2, p3): the angle formed by the three specified points.</li>
 * <li>dihedral(p1, p2, p3, p4): the dihedral angle formed by the four specified
 * points, guaranteed to be in the range [-pi,+pi].</li>
 * </ul>
 *
 * To use this class, create a CustomSummation object, passing the following data to
 * the constructor:
 *
 * <ul>
 * <li>the number of arguments</li>
 * <li>an algebraic expression that defines each term of the sum</li>
 * <li>a map of overall parameter names to default values</li>
 * <li>a list of per-term parameter names</li>
 * <li>the OpenMM::Platform to use for calculations</li>
 * <li>a map of platform-specific property names to values</li>
 * </ul>
 *
 * Then, call addTerm() to define terms of the sum and specify their parameter values.
 * After a term has been added, you can modify its parameters by calling setTerm().
 *
 * As an example, the following code creates a CustomSummation that evaluates a
 * Gaussian mixture in a three-dimensional space. All kernels have the same standard
 * deviation, but different means. A kernel is added for each vertex of a unit cube.
 * Then, the sum is evaluated for a point in the middle of the cube.
 *
 * \verbatim embed:rst:leading-asterisk
 * .. code-block:: cpp
 *
 *     CustomSummation* function = new CustomSummation(
 *         3,
 *         "exp(-((x1-mux)^2+(y1-muy)^2+(z1-muz)^2)/(2*sigma^2))/sqrt(6.2832*sigma^2)",
 *         map<string, double>{{"sigma", 1.0}},
 *         vector<string>{"mux", "muy", "muz"},
 *         Platform::getPlatformByName("CUDA")
 *     );
 *
 *     function->addTerm(vector<double>{0.0, 0.0, 0.0});
 *     function->addTerm(vector<double>{0.0, 0.0, 1.0});
 *     function->addTerm(vector<double>{0.0, 1.0, 0.0});
 *     function->addTerm(vector<double>{0.0, 1.0, 1.0});
 *     function->addTerm(vector<double>{1.0, 0.0, 0.0});
 *     function->addTerm(vector<double>{1.0, 0.0, 1.0});
 *     function->addTerm(vector<double>{1.0, 1.0, 0.0});
 *     function->addTerm(vector<double>{1.0, 1.0, 1.0});
 *
 *     double value = function->evaluate(vector<double>{0.5, 0.5, 0.5});
 *
 * \endverbatim
 *
 * This class also has the ability to compute derivatives of the sum with respect to
 * the arguments.
 *
 * Expressions may involve the operators + (add), - (subtract), * (multiply),
 * / (divide), and ^ (power), and the following functions: sqrt, exp, log, sin, cos,
 * sec, csc, tan, cot, asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max,
 * abs, floor, ceil, step, delta, select.  All trigonometric functions are defined in
 * radians, and log is the natural logarithm.  step(x) = 0 if x is less than 0,
 * 1 otherwise.  delta(x) = 1 if x is 0, 0 otherwise. select(x,y,z) = z if x = 0,
 * y otherwise.
 *
 * This class also supports the functions pointdistance(), pointangle(), and
 * pointdihedral(), which accept 6, 9, and 12 arguments, respectively.  These functions
 * are similar to distance(), angle(), and dihedral(), but their arguments can be any
 * evaluatable expressions rather than the names of predefined points like p1, p2, p3,
 * etc. For example, the following computes the distance from point p1 to the midpoint
 * between p2 and p3.
 *
 * \verbatim embed:rst:leading-asterisk
 * .. code-block:: cpp
 *
 *    string expression = "pointdistance(x1, y1, z1, (x2+x3)/2, (y2+y3)/2, (z2+z3)/2)";
 *
 * \endverbatim
 */

class OPENMM_EXPORT_OPENMM_LAB CustomSummation : public Lepton::CustomFunction {
public:
    /**
     * Construct a new CustomSummation object.
     *
     * @param numArgs                the number of arguments
     * @param expression             the expression for each term in the summation
     * @param overallParameters      the names and default values of the parameters that
     *                               are shared by all terms of the summation. Not to be
     *                               confused with global context parameters
     * @param perTermParameterNames  the names of the parameters that are unique to each
     *                               term of the summation
     * @param platform               the platform that will be used to evaluate the
     *                               summation
     * @param properties             a set of values for platform-specific properties
     */
    CustomSummation(
        int numArgs,
        const std::string &expression,
        const map<string, double> &overallParameters,
        const vector<string> &perTermParameterNames,
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
     * @param arguments    an array of argument values
     */
    double evaluate(const double *arguments) const;
    /**
     * Evaluate the function.
     *
     * @param arguments    a vector of argument values
     */
    double evaluate(const vector<double> &arguments) const;
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
     * Evaluate a derivative of the function.
     *
     * @param arguments    a vector of argument values
     * @param which        the index of the argument for which to evaluate the derivative
     */
    double evaluateDerivative(const vector<double> &arguments, int which) const;
    /**
     * Create a new duplicate of this object on the heap using the "new" operator.
     */
    CustomSummation *clone() const;
    /**
     * Get the expression for each term of the summation.
     *
     * @return         the expression
     */
    const string &getExpression() const { return force->getEnergyFunction(); }
    /**
     * Get the number of overall parameters.
     *
     * @return         the number of overall parameters
     */
    int getNumOverallParameters() const;
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
    /**
     * Get the number of per-term parameters.
     *
     * @return         the number of per-term parameters
     */
    int getNumPerTermParameters() const;
    /**
     * Get the name of a per-term parameter.
     *
     * @param index    the index of the per-term parameter for which to get the name
     * @return         the per-term parameter name
     */
    const std::string &getPerTermParameterName(int index) const;
    /**
     * Get the platform that will be used to evaluate the summation.
     *
     * @return         the platform
     */
    Platform &getPlatform() const { return *platform; }
    /**
     * Get the platform properties.
    */
    const map<string, string> &getPlatformProperties() const;
    /**
     * Add a new term to the summation.
     *
     * @param parameters    the parameters of the term
     * @return              the index of the new term
     */
    int addTerm(vector<double> parameters);
    /**
     * Get the number of terms in the summation.
     *
     * @return         the number of terms
     */
    int getNumTerms() const { return force->getNumBonds(); }
    /**
     * Get the parameters of a term.
     *
     * @param index    the index of the term
     * @return         the parameters of the term
    */
    vector<double> getTerm(int index) const;
    /**
     * Set the parameters of a term.
     *
     * @param index    the index of the term
     * @param parameters    the parameters of the term
     */
    void setTerm(int index, vector<double> parameters);
    /**
     * Get the value of an overall parameter.
     *
     * @param name    the name of the parameter
     * @return        the value of the parameter
     */
    double getParameter(const string &name) const;
    /**
     * Set the value of an overall parameter.
     *
     * @param name    the name of the parameter
     * @param value   the value of the parameter
     */
    void setParameter(const string &name, double value);
private:
    class Evaluator;
    int numArgs;
    vector<int> particles;
    CustomCompoundBondForce *force;
    Evaluator *evaluator;
    Platform *platform;
};

class CustomSummation::Evaluator {
public:
    Evaluator(
        int numArgs,
        CustomCompoundBondForce &force,
        Platform &platform,
        const map<string, string> &properties
    );
    ~Evaluator();
    double evaluate(vector<double> arguments);
    vector<double> evaluateDerivatives(vector<double> arguments);
    void update(CustomCompoundBondForce &force);
    void reset();
    double getParameter(const string &name) const;
    void setParameter(const string &name, double value);
    const map<string, string> &getPlatformProperties() const;
private:
    void setPositions(vector<double> arguments);
    int numArgs;
    Context *context;
    bool contextIsUnchanged;
    vector<Vec3> positions;
    vector<double> latestArguments;
    double value;
    bool valueIsDirty;
    vector<double> derivatives;
    bool derivativesAreDirty;
};

} // namespace OpenMMLab

#endif /*OPENMMLAB_CUSTOMSUMMATION_H_*/
