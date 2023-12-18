import openmmlab as plugin
import numpy as np
import openmm as mm
import pytest
from openmm import unit

cases = [
    ('Reference', ''),
    ('CUDA', 'single'),
    ('CUDA', 'mixed'),
    ('CUDA', 'double'),
    ('OpenCL', 'single'),
    ('OpenCL', 'mixed'),
    ('OpenCL', 'double'),
]

ids = [''.join(case) for case in cases]


def value(x):
    return x/x.unit if unit.is_quantity(x) else x


def ASSERT(cond):
    assert cond


def ASSERT_EQUAL_TOL(expected, found, tol):
    exp = value(expected)
    assert abs(exp - value(found))/max(abs(exp), 1.0) <= tol


def ASSERT_EQUAL_VEC(expected, found, tol):
    ASSERT_EQUAL_TOL(expected.x, found.x, tol)
    ASSERT_EQUAL_TOL(expected.y, found.y, tol)
    ASSERT_EQUAL_TOL(expected.z, found.z, tol)


def assert_forces_and_energy(context, tol):
    state0 = context.getState(getForces=True, getEnergy=True, groups={0})
    state1 = context.getState(getForces=True, getEnergy=True, groups={1})
    for force0, force1 in zip(state0.getForces(), state1.getForces()):
        ASSERT_EQUAL_VEC(force0, force1, tol)
    ASSERT_EQUAL_TOL(state0.getPotentialEnergy(), state1.getPotentialEnergy(), tol)


@pytest.mark.parametrize('platformName, precision', cases, ids=ids)
def testCVs(platformName, precision):
    platform = mm.Platform.getPlatformByName(platformName)
    properties = {} if platformName == 'Reference' else {'Precision': precision}
    system = mm.System()
    system.addParticle(1.0)
    system.addParticle(1.0)
    system.addParticle(1.0)

    # Create a ExtendedCustomCVForce with two collective variables.

    cv = plugin.ExtendedCustomCVForce("v1+3*v2")
    system.addForce(cv)
    v1 = mm.CustomBondForce("2*r")
    v1.addBond(0, 1)
    cv.addCollectiveVariable("v1", v1)
    v2 = mm.CustomBondForce("r")
    v2.addBond(0, 2)
    cv.addCollectiveVariable("v2", v2)
    ASSERT(not cv.usesPeriodicBoundaryConditions())

    # Create a context for evaluating it.

    integrator = mm.VerletIntegrator(1.0)
    context = mm.Context(system, integrator, platform, properties)
    positions = []
    positions.append(mm.Vec3(0, 0, 0))
    positions.append(mm.Vec3(1.5, 0, 0))
    positions.append(mm.Vec3(0, 0, 2.5))
    context.setPositions(positions)

    # Verify the values of the collective variables.

    values = cv.getCollectiveVariableValues(context)
    ASSERT_EQUAL_TOL(2*1.5, values[0], 1e-6)
    ASSERT_EQUAL_TOL(2.5, values[1], 1e-6)

    # Verify the energy and forces.

    state = context.getState(getEnergy=True, getForces=True)
    ASSERT_EQUAL_TOL(2*1.5+3*2.5, state.getPotentialEnergy(), 1e-6)
    ASSERT_EQUAL_VEC(mm.Vec3(2.0, 0.0, 3.0), state.getForces()[0], 1e-6)
    ASSERT_EQUAL_VEC(mm.Vec3(-2.0, 0.0, 0.0), state.getForces()[1], 1e-6)
    ASSERT_EQUAL_VEC(mm.Vec3(0.0, 0.0, -3.0), state.getForces()[2], 1e-6)

@pytest.mark.parametrize('platformName, precision', cases, ids=ids)
def testEnergyParameterDerivatives(platformName, precision):
    platform = mm.Platform.getPlatformByName(platformName)
    properties = {} if platformName == 'Reference' else {'Precision': precision}
    system = mm.System()
    system.addParticle(1.0)
    system.addParticle(1.0)

    # Create a ExtendedCustomCVForce with one collective variable and two global parameters.
    # The CV in turn depends on a global parameter.

    cv = plugin.ExtendedCustomCVForce("v1*g1+g2")
    system.addForce(cv)
    v1 = mm.CustomBondForce("r*g3")
    v1.addGlobalParameter("g3", 2.0)
    v1.addEnergyParameterDerivative("g3")
    v1.addBond(0, 1)
    cv.addCollectiveVariable("v1", v1)
    cv.addGlobalParameter("g1", 1.5)
    cv.addGlobalParameter("g2", -1.0)
    cv.addEnergyParameterDerivative("g2")
    cv.addEnergyParameterDerivative("g1")

    # Create a context for evaluating it.

    integrator = mm.VerletIntegrator(1.0)
    context = mm.Context(system, integrator, platform, properties)
    positions = []
    positions.append(mm.Vec3(0, 0, 0))
    positions.append(mm.Vec3(2.0, 0, 0))
    context.setPositions(positions)

    # Verify the energy and parameter derivatives.

    state = context.getState(getEnergy=True, getParameterDerivatives=True)
    ASSERT_EQUAL_TOL(4.0*1.5-1.0, state.getPotentialEnergy(), 1e-6)
    derivs = state.getEnergyParameterDerivatives()
    ASSERT_EQUAL_TOL(4.0, derivs["g1"], 1e-6)
    ASSERT_EQUAL_TOL(1.0, derivs["g2"], 1e-6)
    ASSERT_EQUAL_TOL(2.0*1.5, derivs["g3"], 1e-6)

@pytest.mark.parametrize('platformName, precision', cases, ids=ids)
def testTabulatedFunction(platformName, precision):
    platform = mm.Platform.getPlatformByName(platformName)
    properties = {} if platformName == 'Reference' else {'Precision': precision}
    xsize = 10
    ysize = 11
    xmin = 0.4
    xmax = 1.1
    ymin = 0.0
    ymax = 0.95
    system = mm.System()
    system.addParticle(1.0)
    integrator = mm.VerletIntegrator(0.01)
    cv = plugin.ExtendedCustomCVForce("fn(x,y)+1")
    v1 = mm.CustomExternalForce("x")
    v1.addParticle(0)
    cv.addCollectiveVariable("x", v1)
    v2 = mm.CustomExternalForce("y")
    v2.addParticle(0)
    cv.addCollectiveVariable("y", v2)
    table = np.empty(xsize*ysize)
    for i in range(xsize):
        for j in range(ysize):
            x = xmin + i*(xmax-xmin)/xsize
            y = ymin + j*(ymax-ymin)/ysize
            table[i+xsize*j] = np.sin(0.25*x)*np.cos(0.33*y)
    tabulated_function = mm.Continuous2DFunction(
        xsize, ysize, table, xmin, xmax, ymin, ymax
    )
    cv.addTabulatedFunction("fn", tabulated_function)
    system.addForce(cv)
    context = mm.Context(system, integrator, platform, properties)
    scale = 1.0
    for i in range(2):
        for x in np.linspace(xmin - 0.15, xmax + 0.15, num=xsize):
            for y in np.linspace(ymin - 0.15, ymax + 0.15, num=ysize):
                positions = [mm.Vec3(x, y, 1.5)]
                context.setPositions(positions)
                state = context.getState(getForces=True, getEnergy=True)
                forces = state.getForces()
                energy = 1
                force = np.zeros(3)
                if (x >= xmin and x <= xmax and y >= ymin and y <= ymax):
                    energy = scale*np.sin(0.25*x)*np.cos(0.33*y)+1
                    force[0] = -scale*0.25*np.cos(0.25*x)*np.cos(0.33*y)
                    force[1] = scale*0.3*np.sin(0.25*x)*np.sin(0.33*y)
                ASSERT_EQUAL_VEC(mm.Vec3(*force), forces[0], 0.1)
                ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), 0.05)

        # Now update the tabulated function, call updateParametersInContext(),
        # and see if it's still correct.

        for j in range(len(table)):
            table[j] *= 2
        tabulated_function.setFunctionParameters(
            xsize, ysize, table, xmin, xmax, ymin, ymax
        )
        cv.updateParametersInContext(context)
        scale *= 2.0

@pytest.mark.parametrize('platformName, precision', cases, ids=ids)
def testReordering(platformName, precision):
    platform = mm.Platform.getPlatformByName(platformName)
    properties = {} if platformName == 'Reference' else {'Precision': precision}

    # Create a larger system with a nonbonded force, since that will trigger atom
    # reordering on the GPU.

    numParticles = 100
    system = mm.System()
    cv = plugin.ExtendedCustomCVForce("2*v2")
    v1 = mm.CustomNonbondedForce("r")
    v1.addPerParticleParameter("a")
    v2 = mm.CustomBondForce("r+1")
    v2.addBond(5, 10)
    cv.addCollectiveVariable("v1", v1)
    cv.addCollectiveVariable("v2", v2)
    cv.setForceGroup(2)
    system.addForce(cv)
    nb = mm.CustomNonbondedForce("r^2")
    nb.addPerParticleParameter("a")
    system.addForce(nb)
    random = np.random.default_rng(0)
    positions = []
    params = np.zeros(1)
    for i in range(numParticles):
        system.addParticle(1.0)
        params[0] = i % 2
        v1.addParticle(params)
        params[0] = 2.0 if i < numParticles/2 else 3.0
        nb.addParticle(params)
        positions.append(mm.Vec3(random.uniform(), random.uniform(), random.uniform())*10)

    # Make sure it works correctly.

    integrator = mm.VerletIntegrator(0.01)
    context = mm.Context(system, integrator, platform, properties)
    context.setPositions(positions)
    state = context.getState(getEnergy=True, getForces=True, groups={2})
    delta = positions[5] - positions[10]
    r = np.linalg.norm(delta)
    ASSERT_EQUAL_TOL(2*(r+1), state.getPotentialEnergy(), 1e-5)
    ASSERT_EQUAL_VEC(-delta*2/r, state.getForces()[5], 1e-5)
    ASSERT_EQUAL_VEC(delta*2/r, state.getForces()[10], 1e-5)
