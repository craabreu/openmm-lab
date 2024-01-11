import openmmlab as plugin
import openmm as mm
import pytest
from openmm import unit

TOL = 1e-5

CASES = [
    ('Reference', ''),
    ('CUDA', 'single'),
    ('CUDA', 'mixed'),
    ('CUDA', 'double'),
    ('OpenCL', 'single'),
    ('OpenCL', 'mixed'),
    ('OpenCL', 'double'),
]

IDS = [''.join(case) for case in CASES]


def value(x):
    return x/x.unit if unit.is_quantity(x) else x


def ASSERT_EQUAL(expected, found, tol=TOL):
    exp = value(expected)
    assert abs(exp - value(found))/max(abs(exp), 1.0) <= tol


@pytest.mark.parametrize('platformName, precision', CASES, ids=IDS)
def testSimpleSummation(platformName, precision):
    platform = mm.Platform.getPlatformByName(platformName)
    properties = {} if platformName == 'Reference' else {'Precision': precision}

    a, b, c, d, e, f, g, h = range(8)
    x1, y1, z1, x2 = range(4)
    overallParameters = {"a": a, "b": b}
    perTermParameters = ("c", "d", "e")

    summation = plugin.CustomSummation(
        4,
        "a*x1+b*y1+c*z1+d*x2+e",
        overallParameters,
        perTermParameters,
        platform,
        properties
    )

    assert summation.getExpression() == "a*x1+b*y1+c*z1+d*x2+e"
    assert summation.getNumArguments() == 4
    assert summation.getOverallParameters() == overallParameters
    assert summation.getPerTermParameters() == perTermParameters
    assert summation.getPlatform().getName() == platformName
    assert summation.getPlatformProperties() == properties

    # Add a term and make sure it is evaluated correctly only after updating.
    summation.addTerm([c, d, e])
    assert summation.getNumTerms() == 1
    args = [x1, y1, z1, x2]
    ASSERT_EQUAL(summation.evaluate(args), 0)
    summation.update()
    ASSERT_EQUAL(summation.evaluate(args), a*x1+b*y1+c*z1+d*x2+e)

    # Add another term and make sure it is evaluated correctly only after updating.
    summation.addTerm([f, g, h])
    assert summation.getNumTerms() == 2
    ASSERT_EQUAL(summation.evaluate(args), a*x1+b*y1+c*z1+d*x2+e)
    summation.update()
    ASSERT_EQUAL(summation.evaluate(args), 2*a*x1+2*b*y1+(c+f)*z1+(d+g)*x2+e+h)

    # Make sure derivatives are evaluated correctly.
    ASSERT_EQUAL(summation.evaluateDerivative(args, 0), 2*a)
    ASSERT_EQUAL(summation.evaluateDerivative(args, 1), 2*b)
    ASSERT_EQUAL(summation.evaluateDerivative(args, 2), c+f)
    ASSERT_EQUAL(summation.evaluateDerivative(args, 3), d+g)

    # Make sure evaluation is correct when new arguments are passed in.
    newArgs = [x1, x1, x1, x1]
    ASSERT_EQUAL(summation.evaluateDerivative(newArgs, 3), d+g)
    ASSERT_EQUAL(summation.evaluateDerivative(newArgs, 2), c+f)
    ASSERT_EQUAL(summation.evaluateDerivative(newArgs, 1), 2*b)
    ASSERT_EQUAL(summation.evaluateDerivative(newArgs, 0), 2*a)
    ASSERT_EQUAL(summation.evaluate(newArgs), (2*(a+b)+c+f+d+g)*x1+e+h)

    # Modify a parameter and make sure evaluation is correct without needing to update.
    anew = 3
    summation.setParameter("a", anew)
    ASSERT_EQUAL(summation.evaluate(newArgs), (2*(anew+b)+c+f+d+g)*x1+e+h)
    ASSERT_EQUAL(summation.evaluateDerivative(args, 0), 2*anew)
