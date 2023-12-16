%module openmmlab

%import(module="openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"
%include <std_string.i>

%{
#include "SlicedNonbondedForce.h"
#include "ExtendedCustomCVForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"

#define SWIG_PYTHON_CAST_MODE
%}

%pythoncode %{
from openmm import unit

__version__ = "@CMAKE_PROJECT_VERSION@"
%}

/*
 * Add units to function outputs.
*/

%pythonappend OpenMMLab::SlicedNonbondedForce::getPMEParametersInContext(
        const OpenMM::Context& context, double& alpha, int& nx, int& ny, int& nz) const %{
    val[0] = unit.Quantity(val[0], 1/unit.nanometers)
%}

%pythonappend OpenMMLab::SlicedNonbondedForce::getLJPMEParametersInContext(
        const OpenMM::Context& context, double& alpha, int& nx, int& ny, int& nz) const %{
    val[0] = unit.Quantity(val[0], 1/unit.nanometers)
%}

/*
 * Convert C++ exceptions to Python exceptions.
*/

%exception {
    try {
        $action
    } catch (std::exception &e) {
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
}

namespace OpenMMLab {

%apply double& OUTPUT {double& alpha};
%apply int& OUTPUT {int& nx};
%apply int& OUTPUT {int& ny};
%apply int& OUTPUT {int& nz};
%apply const std::string& OUTPUT {const std::string& parameter};
%apply int& OUTPUT {int& subset1};
%apply int& OUTPUT {int& subset2};
%apply bool& OUTPUT {bool& includeLJ};
%apply bool& OUTPUT {bool& includeCoulomb};

/**
 * This class implements sliced nonbonded interactions between particles, including a Coulomb force to represent
 * electrostatics and a Lennard-Jones force to represent van der Waals interactions.  It optionally supports
 * periodic boundary conditions and cutoffs for long range interactions.  Lennard-Jones interactions are
 * calculated with the Lorentz-Berthelot combining rule: it uses the arithmetic mean of the sigmas and the
 * geometric mean of the epsilons for the two interacting particles.
 *
 * The particles are classified into :math:`n` disjoint subsets, thus creating :math:`n(n+1)/2` distinct types
 * of particle pairs.  Then, the total potential energy is sliced into this many parts by distinguishing the
 * contributions of such pair types.  Optionally, variables stored as an :OpenMM:`Context` global parameters
 * can multiply the Coulomb and/or Lennard-Jones parts of selected potential energy slices.  Changes in the
 * values of these global parameters via setParameter_ will automatically modify the total potential energy of
 * a SlicedNonbondedForce.  Derivatives with respect to scaling parameters can also be calculated upon request.
 *
 * .. _setParameter: http://docs.openmm.org/latest/api-python/generated/openmm.openmm.Context.html#openmm.openmm.Context.setParameter
 *
 * To use this class, create a SlicedNonbondedForce object, specifying the number :math:`n` of particle subsets.
 * Then, call :func:`addParticle` once for every single particle in the System to define its force field
 * parameters.  Alternatively, you can pass an existing :OpenMM:`NonbondedForce` object when creating a
 * SlicedNonbondedForce object, so that the latter inherits all properties of the former.  After a particle has
 * been added, you can modify its subset index from 0, the default value, to any positive integer lower than
 * :math:`n` by calling :func:`setParticleSubset`.  You can also modify the force field parameters of an added
 * particle by calling :func:`setParticleParameters`.  These two methods will have no effect on Contexts that
 * already exist unless you call :func:`updateParametersInContext`.
 *
 * SlicedNonbondedForce also lets you specify *exceptions* via :func:`addException` and modify their parameters
 * via :func:`setExceptionParameters`.  Exceptions are particular pairs of particles whose interactions should
 * be computed based on different parameters than those defined for the individual particles.  This can be used
 * to exclude certain interactions from the force calculation.
 *
 * By default, all energy slices of a SlicedNonbondedForce are added together, yielding the same total
 * potential energy as a NonbondedForce with equal particle and exception parameters.  In addition, you can
 * select one or more slices and multiply their Coulomb and/or Lennard-Jones contributions by variables
 * referred to as *scaling parameters*.  For this, start by calling :func:`addGlobalParameter` to define
 * a new Context parameter, then pass its name to :func:`addScalingParameter`.  The same scaling parameter
 * can affect multiple slices, and two scaling parameters can affect a single slice, but only if they
 * multiply the Coulomb and Lennard-Jones parts of that slice separately.
 *
 * With :func:`setScalingParameterDerivative`, you can request the derivative with respect to a scaling
 * parameter.  All requested derivatives can then be evaluated for a particular :OpenMM:`State` by calling
 * its getEnergyParameterDerivatives_ method.  For this, such State must have been generated from an
 * :OpenMM:`Context` by passing ``True`` to the keyword ``getParameterDerivatives`` of the Context's
 * getState_ method.
 *
 * .. _getState: http://docs.openmm.org/latest/api-python/generated/openmm.openmm.Context.html#openmm.openmm.Context.getState
 *
 * .. _getEnergyParameterDerivatives: http://docs.openmm.org/latest/api-python/generated/openmm.openmm.State.html#openmm.openmm.State.getEnergyParameterDerivatives
 *
 * Many molecular force fields omit Coulomb and Lennard-Jones interactions between particles separated by one
 * or two bonds, while using modified parameters for those separated by three bonds (known as "1-4 interactions").
 * This class provides a convenience method for this case called :func:`createExceptionsFromBonds`.  You pass to it
 * a list of bonds and the scale factors to use for 1-4 interactions.  It identifies all pairs of particles which
 * are separated by 1, 2, or 3 bonds, then automatically creates exceptions for them.
 *
 * When using a cutoff, by default Lennard-Jones interactions are sharply truncated at the cutoff distance.
 * Optionally you can instead use a switching function to make the interaction smoothly go to zero over a finite
 * distance range.  To enable this, call :func:`setUseSwitchingFunction`.  You must also call :func:`setSwitchingDistance`
 * to specify the distance at which the interaction should begin to decrease.  The switching distance must be
 * less than the cutoff distance.
 *
 * Another optional feature of this class (enabled by default) is to add a contribution to the energy which approximates
 * the effect of all Lennard-Jones interactions beyond the cutoff in a periodic system.  When running a simulation
 * at constant pressure, this can improve the quality of the result.  Call :func:`setUseDispersionCorrection` to set whether
 * this should be used.
 *
 * In some applications, it is useful to be able to inexpensively change the parameters of small groups of particles.
 * Usually this is done to interpolate between two sets of parameters.  For example, a titratable group might have
 * two states it can exist in, each described by a different set of parameters for the atoms that make up the
 * group.  You might then want to smoothly interpolate between the two states.  This is done by first calling
 * :func:`addGlobalParameter` to define a Context parameter, then :func:`addParticleParameterOffset` to create a "parameter offset"
 * that depends on the Context parameter.  Each offset defines the following:
 *
 * * A Context parameter used to interpolate between the states.
 * * A single particle whose parameters are influenced by the Context parameter.
 * * Three scale factors (chargeScale, sigmaScale, and epsilonScale) that specify how the Context parameter affects the particle.
 *
 * The "effective" parameters for a particle (those used to compute forces) are given by
 *
 * .. code-block:: python
 *
 *    charge = baseCharge + param*chargeScale
 *    sigma = baseSigma + param*sigmaScale
 *    epsilon = baseEpsilon + param*epsilonScale
 *
 * where the ``base`` values are the ones specified by :func:`addParticle` and ``param`` is the current value
 * of the Context parameter.  A single Context parameter can apply offsets to multiple particles,
 * and multiple parameters can be used to apply offsets to the same particle.  Parameters can also be used
 * to modify exceptions in exactly the same way by calling :func:`addExceptionParameterOffset`.  A context
 * parameter cannot be used simultaneously to apply offsets and scale energy terms.
 */
class SlicedNonbondedForce : public OpenMM::NonbondedForce {
public:
    /**
     * Create a SlicedNonbondedForce.
     *
     * Parameters
     * ----------
     *     numSubsets : int
     *         the number of particle subsets
     */
    SlicedNonbondedForce(int numSubsets);
    /**
     * Create a SlicedNonbondedForce having the properties of an existing :OpenMM:`NonbondedForce`.
     *
     * Parameters
     * ----------
     *     force : :OpenMM:`NonbondedForce`
     *         the NonbondedForce object from which to instantiate this SlicedNonbondedForce
     *     numSubsets : int
     *         the number of particle subsets
     */
    SlicedNonbondedForce(const OpenMM::NonbondedForce& force, int numSubsets);
    /**
     * Get the parameters being used for PME in a particular Context.  Because some platforms have restrictions
     * on the allowed grid sizes, the values that are actually used may be slightly different from those
     * specified with setPMEParameters(), or the standard values calculated based on the Ewald error tolerance.
     * See the manual for details.
     *
     * Parameters
     * ----------
     *     context : Context
     *         the Context for which to get the parameters
     *
     * Returns
     * -------
     *     alpha : double
     *         the separation parameter, measured in :math:`nm^{-1}`
     *     nx : int
     *         the number of grid points along the X axis
     *     ny : int
     *         the number of grid points along the Y axis
     *     nz : int
     *         the number of grid points along the Z axis
     */
    void getPMEParametersInContext(const OpenMM::Context& context, double& alpha, int& nx, int& ny, int& nz) const;
    /**
     * Get the PME parameters being used for the dispersion term for LJPME in a particular Context.  Because some
     * platforms have restrictions on the allowed grid sizes, the values that are actually used may be slightly different
     * from those specified with setPMEParameters(), or the standard values calculated based on the Ewald error tolerance.
     * See the manual for details.
     *
     * Parameters
     * ----------
     *     context : Context
     *         the Context for which to get the parameters
     *
     * Returns
     * -------
     *     alpha : double
     *         the separation parameter, measured in :math:`nm^{-1}`
     *     nx : int
     *         the number of grid points along the X axis
     *     ny : int
     *         the number of grid points along the Y axis
     *     nz : int
     *         the number of grid points along the Z axis
     */
    void getLJPMEParametersInContext(const OpenMM::Context& context, double& alpha, int& nx, int& ny, int& nz) const;
    /**
     * Update the particle and exception parameters in a Context to match those stored in this Force object.  This method
     * provides an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setParticleParameters() and setExceptionParameters() to modify this object's parameters, then call
     * updateParametersInContext() to copy them over to the Context.
     *
     * This method has several limitations.  The only information it updates is the parameters of particles and exceptions.
     * All other aspects of the Force (the nonbonded method, the cutoff distance, etc.) are unaffected and can only be
     * changed by reinitializing the Context.  Furthermore, only the chargeProd, sigma, and epsilon values of an exception
     * can be changed; the pair of particles involved in the exception cannot change.  Finally, this method cannot be used
     * to add new particles or exceptions, only to change the parameters of existing ones.
     *
     * Parameters
     * ----------
     *     context : Context
     *         the Context in which to update the parameters
     */
    void updateParametersInContext(OpenMM::Context& context);
    /**
     * Get the name of the method used for handling long range nonbonded interactions.
     */
    std::string getNonbondedMethodName() const;
    /**
     * Get the specified number of particle subsets.
     */
    int getNumSubsets() const;
    /**
     * Get the number of slices determined by the specified number of particle subsets.
     */
    int getNumSlices() const;
    /**
     * Get the subset to which a particle belongs.
     *
     * Parameters
     * ----------
     *     index : int
     *         the index of the particle for which to get the subset
     */
    int getParticleSubset(int index) const;
    /**
     * Set the subset of a particle.
     *
     * Parameters
     * ----------
     *     index : int
     *         the index of the particle for which to set the subset
     *     subset : int
     *         the subset to which this particle belongs
     */
    void setParticleSubset(int index, int subset);
  	/**
     * Add a scaling parameter to multiply a particular Coulomb slice. Its value will scale the
     * Coulomb interactions between particles of a subset 1 with those of another (or the same)
     * subset 2. The order of subset definition is irrelevant.
     *
     * Parameters
     * ----------
     *     parameter : str
     *         the name of the global parameter.  It must have already been added
     *         with :func:`addGlobalParameter`. Its value can be modified at any time by
     *         calling `setParameter()` on the :OpenMM:`Context`
     *     subset1 : int
     *         the index of a particle subset.  Legal values are between 0 and the result of
     *         :func:`getNumSubsets`
     *     subset2 : int
     *         the index of a particle subset.  Legal values are between 0 and the result of
     *         :func:`getNumSubsets`
     *     includeCoulomb : bool
     *         whether this scaling parameter applies to Coulomb interactions
     *     includeLJ : bool
     *         whether this scaling parameter applies to Lennard-Jones interactions
     *
     * Returns
     * -------
     *     index : int
     *         the index of scaling parameter that was added
     */
    int addScalingParameter(const std::string& parameter, int subset1, int subset2, bool includeCoulomb, bool includeLJ);
    /**
     * Get the number of scaling parameters.
     */
    int getNumScalingParameters() const;
  	/**
     * Get the scaling parameter applied to a particular nonbonded slice.
     *
     * Parameters
     * ----------
     *     index : int
     *         the index of the scaling parameter to query, as returned by :func:`addScalingParameter`
     *
     * Returns
     * -------
     *     parameter : str
     *         the name of the global parameter
     *     subset1 : int
     *         the smallest index of the two particle subsets
     *     subset2 : int
     *         the largest index of the two particle subsets
     *     includeCoulomb : bool
     *         whether this scaling parameter applies to Coulomb interactions
     *     includeLJ : bool
     *         whether this scaling parameter applies to Lennard-Jones interactions
     */
    void getScalingParameter(int index, std::string& parameter, int& subset1, int& , bool& includeCoulomb, bool& includeLJ) const;
 	/**
     * Modify an added scaling parameter.
     *
     * Parameters
     * ----------
     *     index : int
     *         the index of the scaling parameter to modify, as returned by
     *         :func:`addExceptionChargeOffset`
     *     parameter : str
     *         the name of the global parameter.  It must have already been added
     *         with :func:`addGlobalParameter`. Its value can be modified at any time by
     *         calling `setParameter()` on the :OpenMM:`Context`
     *     subset1 : int
     *         the index of a particle subset.  Legal values are between 0 and the result of
     *         :func:`getNumSubsets`
     *     subset2 : int
     *         the index of a particle subset.  Legal values are between 0 and the result of
     *         :func:`getNumSubsets`
     *     includeCoulomb : bool
     *         whether this scaling parameter applies to Coulomb interactions
     *     includeLJ : bool
     *         whether this scaling parameter applies to Lennard-Jones interactions
     */
    void setScalingParameter(int index, const std::string& parameter, int subset1, int subset2, bool includeCoulomb, bool includeLJ);
    /**
     * Request the derivative of this Force's energy with respect to a scaling parameter. This
     * can be used to obtain the sum of particular energy slices. The parameter must have already
     * been added with :func:`addGlobalParameter` and :func:`addScalingParameter`.
     *
     * Parameters
     * ----------
     *     parameter : str
     *         the name of the parameter
     *
     * Returns
     * -------
     *     index : int
     *         the index of scaling parameter derivative that was added
     */
    int addScalingParameterDerivative(const std::string& parameter);
    /**
     * Get the number of requested scaling parameter derivatives.
     */
    int getNumScalingParameterDerivatives() const;
    /**
     * Get the name of the global parameter associated with a requested scaling parameter
     * derivative.
     *
     * Parameters
     * ----------
     *     index : int
     *         the index of the parameter derivative, between 0 and the result of
     *         :func:`getNumScalingParameterDerivatives`
     */
    const std::string& getScalingParameterDerivativeName(int index) const;
    /**
     * Set the name of the global parameter to associate with a requested scaling parameter
     * derivative.
     *
     * Parameters
     * ----------
     *     index : int
     *         the index of the parameter derivative, between 0 and getNumScalingParameterDerivatives`
     *     parameter : str
     *         the name of the parameter
     */
    void setScalingParameterDerivative(int index, const std::string& parameter);
	/**
     * Get whether to use CUDA Toolkit's cuFFT library when executing in the CUDA platform.
     * The default value is `False`.
     */
    bool getUseCudaFFT() const;
 	/**
     * Set whether whether to use CUDA Toolkit's cuFFT library when executing in the CUDA platform.
     * This choice has no effect when using other platforms or when the CUDA Toolkit is version 7.0
     * or older.
     *
     * Parameters
     * ----------
     *     use : bool
     *         whether to use the cuFFT library
     */
    void setUseCuFFT(bool use);

    /*
     * Add methods for casting a Force to a SlicedNonbondedForce.
    */

    %extend {
        static OpenMMLab::SlicedNonbondedForce& cast(OpenMM::Force& force) {
            return dynamic_cast<OpenMMLab::SlicedNonbondedForce&>(force);
        }

        static bool isinstance(OpenMM::Force& force) {
            return (dynamic_cast<OpenMMLab::SlicedNonbondedForce*>(&force) != NULL);
        }
    }
};

%clear double& alpha;
%clear int& nx;
%clear int& ny;
%clear int& nz;
%clear std::string& parameter;
%clear int& subset1;
%clear int& subset2;
%clear bool& includeLJ;
%clear bool& includeCoulomb;

/**
 * This class supports energy functions that depend on collective variables.  To use it,
 * you define a set of collective variables (scalar valued functions that depend on the
 * particle positions), and an algebraic expression for the energy as a function of the
 * collective variables.  The expression also may involve tabulated functions, and may
 * depend on arbitrary global parameters.
 *
 * Each collective variable is defined by a Force object.  The Force's potential energy
 * is computed, and that becomes the value of the variable.  This provides enormous
 * flexibility in defining collective variables, especially by using custom forces.
 * Anything that can be computed as a potential function can also be used as a collective
 * variable.
 *
 * To use this class, create a ExtendedCustomCVForce object, passing an algebraic expression to the
 * constructor that defines the potential energy.  Then call addCollectiveVariable() to define
 * collective variables and addGlobalParameter() to define global parameters.  The values
 * of global parameters may be modified during a simulation by calling Context::setParameter().
 *
 * This class also has the ability to compute derivatives of the potential energy with respect to global parameters.
 * Call addEnergyParameterDerivative() to request that the derivative with respect to a particular parameter be
 * computed.  You can then query its value in a Context by calling getState() on it.
 *
 * Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and ^ (power), and the following
 * functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max, abs, floor, ceil, step, delta, select.  All trigonometric functions
 * are defined in radians, and log is the natural logarithm.  step(x) = 0 if x is less than 0, 1 otherwise.  delta(x) = 1 if x is 0, 0 otherwise.
 * select(x,y,z) = z if x = 0, y otherwise.
 *
 * In addition, you can call addTabulatedFunction() to define a new function based on tabulated values.  You specify the function by
 * creating a TabulatedFunction object.  That function can then appear in the expression.
 */
class ExtendedCustomCVForce : public OpenMM::Force {
public:
    /**
     * Create a ExtendedCustomCVForce.
     *
     * Parameters
     * ----------
     * energy : str
     *     an algebraic expression giving the energy of the system as a function
     *     of the collective variables and global parameters
     */
    explicit ExtendedCustomCVForce(const std::string& energy);
    /**
     * Get the number of collective variables that the interaction depends on.
     */
    int getNumCollectiveVariables() const;
    /**
     * Get the number of global parameters that the interaction depends on.
     */
    int getNumGlobalParameters() const;
    /**
     * Get the number of global parameters with respect to which the derivative of the energy
     * should be computed.
     */
    int getNumEnergyParameterDerivatives() const;
    /**
     * Get the number of tabulated functions that have been defined.
     */
    int getNumTabulatedFunctions() const;
    /**
     * Get the algebraic expression that gives the energy of the system
     */
    const std::string& getEnergyFunction() const;
    /**
     * Set the algebraic expression that gives the energy of the system
     */
    void setEnergyFunction(const std::string& energy);
    /**
     * Add a collective variable that the force may depend on.  The collective variable
     * is represented by a Force object, which should have been created on the heap with the
     * "new" operator.  The ExtendedCustomCVForce takes over ownership of it, and deletes the Force when the
     * ExtendedCustomCVForce itself is deleted.
     *
     * Parameters
     * ----------
     * name : str
     *     the name of the collective variable, as it will appear in the energy expression
     * variable : Force
     *     the collective variable, represented by a Force object.  The value of the
     *     variable is the energy computed by the Force.
     *
     * Returns
     * -------
     * int
     *     the index within the Force of the variable that was added
     */
    int addCollectiveVariable(const std::string& name, Force* variable);
    /**
     * Get the name of a collective variable.
     *
     * Parameters
     * ----------
     * index : int
     *     the index of the collective variable for which to get the name
     *
     * Returns
     * -------
     * str
     *     the variable name
     */
    const std::string& getCollectiveVariableName(int index) const;
    /**
     * Get a writable reference to the Force object that computes a collective variable.
     *
     * Parameters
     * ----------
     * index : int
     *     the index of the collective variable to get
     *
     * Returns
     * -------
     * Force
     *     the Force object
     */
    Force& getCollectiveVariable(int index);
    /**
     * Get a const reference to the Force object that computes a collective variable.
     *
     * Parameters
     * ----------
     * index : int
     *     the index of the collective variable to get
     *
     * Returns
     * -------
     * const Force
     *     the Force object
     */
    const Force& getCollectiveVariable(int index) const;
    /**
     * Add a new global parameter that the interaction may depend on.  The default value provided to
     * this method is the initial value of the parameter in newly created Contexts.  You can change
     * the value at any time by calling setParameter() on the Context.
     *
     * Parameters
     * ----------
     * name : str
     *     the name of the parameter
     * default_value : float
     *     the default value of the parameter
     *
     * Returns
     * -------
     * int
     *     the index of the parameter that was added
     */
    int addGlobalParameter(const std::string& name, double defaultValue);
    /**
     * Get the name of a global parameter.
     *
     * Parameters
     * ----------
     * index : int
     *     the index of the parameter for which to get the name
     *
     * Returns
     * -------
     * str
     *     the parameter name
     */
    const std::string& getGlobalParameterName(int index) const;
    /**
     * Set the name of a global parameter.
     *
     * Parameters
     * ----------
     * index : int
     *     the index of the parameter for which to set the name
     * name : str
     *     the name of the parameter
     */
    void setGlobalParameterName(int index, const std::string& name);
    /**
     * Get the default value of a global parameter.
     *
     * Parameters
     * ----------
     * index : int
     *     the index of the parameter for which to get the default value
     *
     * Returns
     * -------
     * float
     *     the parameter default value
     */
    double getGlobalParameterDefaultValue(int index) const;
    /**
     * Set the default value of a global parameter.
     *
     * Parameters
     * ----------
     * index : int
     *     the index of the parameter for which to set the default value
     * default_value : float
     *     the default value of the parameter
     *
     * Returns
     * -------
     * float
     *     the parameter default value
     */
    void setGlobalParameterDefaultValue(int index, double defaultValue);
    /**
     * Request that this Force compute the derivative of its energy with respect to a global parameter.
     * The parameter must have already been added with addGlobalParameter().
     *
     * Parameters
     * ----------
     * name : str
     *     the name of the parameter
     */
    void addEnergyParameterDerivative(const std::string& name);
    /**
     * Get the name of a global parameter with respect to which this Force should compute the
     * derivative of the energy.
     *
     * Parameters
     * ----------
     * index : int
     *     the index of the parameter derivative, between 0 and getNumEnergyParameterDerivatives()
     *
     * Returns
     * -------
     * str
     *     the parameter name
     */
    const std::string& getEnergyParameterDerivativeName(int index) const;
    /**
     * Add a tabulated function that may appear in the energy expression.
     *
     * Parameters
     * ----------
     * name : str
     *     the name of the function
     * function : TabulatedFunction
     *     a TabulatedFunction object defining the function
     *
     * Returns
     * -------
     * int
     *     the index of the function that was added
     */
    int addTabulatedFunction(const std::string& name, TabulatedFunction* function);
    /**
     * Get a const reference to a tabulated function that may appear in the energy expression.
     *
     * Parameters
     * ----------
     * index : int
     *     the index of the function to get
     *
     * Returns
     * -------
     * const TabulatedFunction
     *     the TabulatedFunction object defining the function
     */
    const TabulatedFunction& getTabulatedFunction(int index) const;
    /**
     * Get a reference to a tabulated function that may appear in the energy expression.
     *
     * Parameters
     * ----------
     * index : int
     *     the index of the function to get
     *
     * Returns
     * -------
     * TabulatedFunction
     *     the TabulatedFunction object defining the function
     */
    TabulatedFunction& getTabulatedFunction(int index);
    /**
     * Get the name of a tabulated function that may appear in the energy expression.
     *
     * Parameters
     * ----------
     * index : int
     *     the index of the function to get
     *
     * Returns
     * -------
     * str
     *     the name of the function as it appears in expressions
     */
    const std::string& getTabulatedFunctionName(int index) const;
    /**
     * Get the current values of the collective variables in a Context.
     *
     * Parameters
     * ----------
     * context : Context
     *     the Context
     *
     * Returns
     * -------
     * std::vector<double>
     *     the values of the collective variables
     */
#if OPENMM_VERSION_MAJOR >= 8
    void getCollectiveVariableValues(OpenMM::Context& context, std::vector<double>& values) const;
#else
    void getCollectiveVariableValues(OpenMM::Context& context, std::vector<double>& values);
#endif
    /**
     * Get the inner Context used for evaluating collective variables.
     *
     * When you create a Context for a System that contains a ExtendedCustomCVForce, internally
     * it creates a new System, adds the Forces that define the CVs to it, creates a new
     * Context for that System, and uses it to evaluate the variables.  In most cases you
     * can ignore all of this.  It is just an implementation detail.  However, there are
     * a few cases where you need to directly access that internal Context.  For example,
     * if you want to modify one of the Forces that defines a collective variable and
     * call updateParametersInContext() on it, you need to pass that inner Context to it.
     * This method returns a reference to it.
     *
     * Parameters
     * ----------
     * context : Context
     *     the Context
     *
     * Returns
     * -------
     * Context
     *     the inner Context used to evaluate the collective variables
     */
    OpenMM::Context& getInnerContext(OpenMM::Context& context);
    /**
     * Update the tabulated function parameters in a Context to match those stored in this Force object.  This method
     * provides an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call getTabulatedFunction(index).setFunctionParameters() to modify this object's parameters, then call
     * updateParametersInContext() to copy them over to the Context.
     *
     * This method is very limited.  The only information it updates is the parameters of tabulated functions.
     * All other aspects of the Force (the energy expression, the set of collective variables, etc.) are unaffected and can
     * only be changed by reinitializing the Context.
     */
    void updateParametersInContext(OpenMM::Context& context);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * Returns
     * -------
     * bool
     *     true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const;

    /*
     * Add methods for casting a Force to an ExtendedCustomCVForce.
    */

    %extend {
        static OpenMMLab::ExtendedCustomCVForce& cast(OpenMM::Force& force) {
            return dynamic_cast<OpenMMLab::ExtendedCustomCVForce&>(force);
        }

        static bool isinstance(OpenMM::Force& force) {
            return (dynamic_cast<OpenMMLab::ExtendedCustomCVForce*>(&force) != NULL);
        }
    }
};

}
