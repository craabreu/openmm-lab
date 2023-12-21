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
#include "openmm/internal/AssertionUtilities.h"
#include <iostream>
#include <string>

using namespace OpenMMLab;
using namespace OpenMM;
using namespace std;


void testSimpleSummation() {
    const double a = 1.0, b = 2.0, c = 3.0, d = 4.0, e = 5.0, f = 6.0, g = 7.0, h = 8.0;
    const double x1 = 1.0, y1 = 2.0, z1 = 3.0, x2 = 4.0;
    map<string, double> overallParameters = {{"a", a}, {"b", b}};
    vector<string> perTermParameters = {"c", "d", "e"};
    CustomSummation summation(
        4,
        "a*x1+b*y1+c*z1+d*x2+e",
        overallParameters,
        perTermParameters,
        platform,
        properties
    );
    ASSERT_EQUAL(summation.getExpression(), "a*x1+b*y1+c*z1+d*x2+e");
    ASSERT_EQUAL(summation.getNumArguments(), 4);
    ASSERT_EQUAL(summation.getNumPerTermParameters(), 3);
    ASSERT_EQUAL(summation.getNumTerms(), 0);
    ASSERT_EQUAL(summation.getNumOverallParameters(), 2);
    ASSERT_EQUAL(summation.getOverallParameterName(0), "a");
    ASSERT_EQUAL(summation.getOverallParameterDefaultValue(0), a);
    ASSERT_EQUAL(summation.getOverallParameterName(1), "b");
    ASSERT_EQUAL(summation.getOverallParameterDefaultValue(1), b);

    summation.addTerm(vector<double>{c, d, e});
    ASSERT_EQUAL(summation.getNumTerms(), 1);
    double *args = new double[4]{x1, y1, z1, x2};
    ASSERT_EQUAL(summation.evaluate(args), a*x1+b*y1+c*z1+d*x2+e);
    summation.addTerm(vector<double>{f, g, h});
    ASSERT_EQUAL(summation.getNumTerms(), 2);
    int *derivOrder = new int[4]{0, 0, 0, 0};
    ASSERT_EQUAL(summation.evaluate(args), 2*a*x1+2*b*y1+(c+f)*z1+(d+g)*x2+e+h);
    derivOrder[0] = 1;
    ASSERT_EQUAL(summation.evaluateDerivative(args, derivOrder), 2*a);
    derivOrder[0] = 0;
    derivOrder[1] = 1;
    ASSERT_EQUAL(summation.evaluateDerivative(args, derivOrder), 2*b);
    derivOrder[1] = 0;
    derivOrder[2] = 1;
    ASSERT_EQUAL(summation.evaluateDerivative(args, derivOrder), c+f);
    derivOrder[2] = 0;
    derivOrder[3] = 1;
    ASSERT_EQUAL(summation.evaluateDerivative(args, derivOrder), d+g);

    vector<double> newArgs({x1, x1, x1, x1});
    ASSERT_EQUAL(summation.evaluateDerivative(newArgs, 3), d+g);
    ASSERT_EQUAL(summation.evaluateDerivative(newArgs, 2), c+f);
    ASSERT_EQUAL(summation.evaluateDerivative(newArgs, 1), 2*b);
    ASSERT_EQUAL(summation.evaluateDerivative(newArgs, 0), 2*a);
    ASSERT_EQUAL(summation.evaluate(newArgs), (2*(a+b)+c+f+d+g)*x1+e+h);

    delete[] args;
    delete[] derivOrder;
}

void testCloning() {
    const double a = 1.0, b = 2.0, c = 3.0, d = 4.0, e = 5.0, f = 6.0, g = 7.0;
    const double x1 = 1.0, y1 = 2.0;
    map<string, double> overallParameters = {{"a", a}};
    vector<string> perTermParameters = {"b", "c"};
    CustomSummation summation(
        2,
        "a*x1+b*y1+c",
        overallParameters,
        perTermParameters,
        platform,
        properties
    );
    summation.addTerm(vector<double>{b, c});
    summation.addTerm(vector<double>{d, e});
    summation.addTerm(vector<double>{f, g});
    CustomSummation *copy = summation.clone();
    ASSERT_EQUAL(copy->getExpression(), "a*x1+b*y1+c");
    ASSERT_EQUAL(copy->getNumArguments(), 2);
    ASSERT_EQUAL(copy->getNumPerTermParameters(), 2);
    ASSERT_EQUAL(copy->getNumTerms(), 3);
    ASSERT_EQUAL(copy->getNumOverallParameters(), 1);
    ASSERT_EQUAL(copy->getOverallParameterName(0), "a");
    ASSERT_EQUAL(copy->getOverallParameterDefaultValue(0), a);
    ASSERT_EQUAL(
        copy->evaluate(vector<double>{x1, y1}),
        3*a*x1+(b+d+f)*y1+c+e+g
    )

    delete copy;
}

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testSimpleSummation();
        testCloning();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
