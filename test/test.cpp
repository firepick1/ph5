#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "FireLog.h"
#include "FireUtils.hpp"
#include "ph5.hpp"
#include "version.h"

using namespace std;
using namespace ph5;

void test_Complex() {
    cout << "TEST	: test_Complex()" << endl;

    Complex<float> c1(1, 2);
    ASSERTEQUAL(c1.Re(), 1);
    ASSERTEQUAL(c1.Im(), 2);

    Complex<float> c5(3, 4);
    ASSERTEQUAL(5, c5.modulus());

    Complex<float> c15 = c1 + c5;
    ASSERT(c15 == Complex<float>(4, 6));
    ASSERT(c1 + 1 == Complex<float>(2, 2));

    Complex<float> c13(1, 3);
    Complex<float> c24(2, 4);
    Complex<float> c1324 = c13 - c24;
    ASSERT((c13 - c24) == Complex<float>(-1, -1));
    ASSERT((c24 - c13) == Complex<float>(1, 1));
    ASSERT((c13 - 1) == Complex<float>(0, 3));

    ASSERT((c13 * c24) == Complex<float>(1 * 2 - 3 * 4, 1 * 4 + 3 * 2));
    ASSERT((c24 * c13) == Complex<float>(1 * 2 - 3 * 4, 1 * 4 + 3 * 2));
    ASSERT((c13 * 2) == Complex<float>(2, 6));

    ASSERT(c13.conj() == Complex<float>(1, -3));

    Complex<float> cr = c13.recip();
    float epsilon = 0.000001;
    ASSERTEQUALT(cr.Re(), 0.1, epsilon);
    ASSERTEQUALT(cr.Im(), -0.3, epsilon);
    Complex<float> crr = cr.recip();
    ASSERTEQUALT(crr.Re(), 1, epsilon);
    ASSERTEQUALT(crr.Im(), 3, epsilon);

    Complex<float> c13div24 = c13 / c24;
    ASSERTEQUALT(c13div24.Re(), 0.7, epsilon);
    ASSERTEQUALT(c13div24.Im(), 0.1, epsilon);

    Complex<float> sum1;
    sum1.add(c13);
    ASSERT(sum1 == Complex<float>(1, 3));
    sum1.add(c24);
    ASSERT(sum1 == Complex<float>(3, 7));

    ASSERT(Complex<float>(25).sqrt() == Complex<float>(5));
    ASSERT(Complex<float>(-1).sqrt() == Complex<float>(0, 1));
    ASSERT(Complex<float>(3, 4).sqrt() == Complex<float>(2, 1));

    ASSERTEQUALS("0", Complex<float>().stringify().c_str());
    ASSERTEQUALS("1", Complex<float>(1).stringify().c_str());
    ASSERTEQUALS("i", Complex<float>(0, 1).stringify().c_str());
    ASSERTEQUALS("-i", Complex<float>(0, -1).stringify().c_str());
    ASSERTEQUALS("2i", Complex<float>(0, 2).stringify().c_str());
    ASSERTEQUALS("-2i", Complex<float>(0, -2).stringify().c_str());
    ASSERTEQUALS("1+i", Complex<float>(1, 1).stringify().c_str());
    ASSERTEQUALS("1-i", Complex<float>(1, -1).stringify().c_str());
    ASSERTEQUALS("1+2i", Complex<float>(1, 2).stringify().c_str());
    ASSERTEQUALS("1-2i", Complex<float>(1, -2).stringify().c_str());
    ASSERTEQUALS("1.2-5.7i", Complex<float>(1.2345, -5.6789).stringify(1).c_str());
    ASSERTEQUALS("1.0+1.0i", Complex<float>(0.99, 0.99).stringify(1).c_str());
    ASSERTEQUALS("1.0-1.0i", Complex<float>(0.99, -0.99).stringify(1).c_str());
    ASSERTEQUALS("0.99-0.99i", Complex<float>(0.99, -0.99).stringify(2).c_str());

    c13.assertEqualT(Complex<float>(1.0009, 2.9991), 0.001);
    cout << "TEST	:   test_Complex() OK" << endl;
}

void test_Bernstein() {
    float  epsilon = 0.00001;
    cout << "TEST	: test_Bernstein()" << endl;
    ASSERTEQUALT(0, Bernstein5<float>(5, 0), epsilon);
    ASSERTEQUALT(0.03125, Bernstein5<float>(5, 0.5), epsilon);
    ASSERTEQUALT(0.15625, Bernstein5<float>(1, 0.5), epsilon);
    cout << "TEST	:   test_Bernstein() OK" << endl;
}

PH5Curve<float> ph_arc() {
    vector<Complex<float> > q;
    vector<Complex<float> > z;
    q.push_back(Complex<float>(-1, 1));
    q.push_back(Complex<float>(0, 2));
    q.push_back(Complex<float>(1, 1));
    z.push_back(Complex<float>());
    z.push_back(Complex<float>(1.124171968973597, 0.444771808762066));
    z.push_back(Complex<float>(1.124171968973597, -0.444771808762066));
    return PH5Curve<float>(z, q);
}

void test_PH5Curve() {
    cout << "TEST	: test_PH5Curve()" << endl;
    PH5Curve<float> ph(ph_arc());

    long msStart = millis();
    float epsilon = 0.00001;
    int ITER = 10000;
    for (int i = 0; i < ITER; i++) {
        Complex<float> c;
        c = ph.r(0);
        c.assertEqualT(Complex<float>(-1, 1), epsilon);
        c = ph.r(0.25);
        c.assertEqualT(Complex<float>(-0.598911, 1.75), epsilon);
        c = ph.r(0.5);
        c.assertEqualT(Complex<float>(0, 2), epsilon);
        c = ph.r(0.75);
        c.assertEqualT(Complex<float>(0.598911, 1.75), epsilon);
        c = ph.r(1.0);
        c.assertEqualT(Complex<float>(1, 1), epsilon);
    }
    long msElapsed = millis() - msStart;
    cout << "TEST	:   average execution time for r(): " << msElapsed / (ITER / 5000) << "us" << endl;

    ASSERTEQUALT(0.00000, ph.s(0.0), epsilon);
    ASSERTEQUALT(1.52753, ph.s(0.5), epsilon);
    ASSERTEQUALT(3.05505, ph.s(1.0), epsilon);

    ASSERTEQUALT(4.11010, ph.sigma(0.0), epsilon);
    ASSERTEQUALT(2.78074, ph.sigma(0.3), epsilon);
    ASSERTEQUALT(2.52753, ph.sigma(0.5), epsilon);
    ASSERTEQUALT(2.78074, ph.sigma(0.7), epsilon);
    ASSERTEQUALT(4.11010, ph.sigma(1.0), epsilon);

    epsilon = 0.001;
    ph.rprime(0.00).assertEqualT(Complex<float>(0.945, 4.000), epsilon);
    ph.rprime(0.25).assertEqualT(Complex<float>(2.132, 2.000), epsilon);
    ph.rprime(0.50).assertEqualT(Complex<float>(2.528, 0.000), epsilon);
    ph.rprime(0.75).assertEqualT(Complex<float>(2.132, -2.000), epsilon);
    ph.rprime(1.00).assertEqualT(Complex<float>(0.945, -4.000), epsilon);

    cout << "TEST	:   test_PH5Curve() OK " << endl;
}

void test_PHFeed() {
    cout << "TEST	: test_PH5Feed()" << endl;
    PH5Curve<float> ph(ph_arc());
    PHFeed<float> phf(ph, 100, 0.01, 0, 100, 0);
    float E = 0;
    float epsilon = 0.001;
    ASSERTEQUALT(0, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<float>(-1.000, 1.000), epsilon);
    E = phf.Ekt(E, 0.1);
    ASSERTEQUALT(0.039, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<float>(-0.990, 1.038), 0.001);
    E = phf.Ekt(E, 0.2);
    ASSERTEQUALT(0.313, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<float>(-0.904, 1.298), 0.001);
    E = phf.Ekt(E, 0.3);
    ASSERTEQUALT(0.716, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<float>(-0.699, 1.643), 0.001);
    E = phf.Ekt(E, 0.4);
    ASSERTEQUALT(1.122, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<float>(-0.389, 1.901), 0.001);
    E = phf.Ekt(E, 0.5);
    ASSERTEQUALT(1.527, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<float>(-0.000, 2.000), 0.001);
    E = phf.Ekt(E, 0.6);
    ASSERTEQUALT(1.933, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<float>(0.389, 1.901), 0.001);
    E = phf.Ekt(E, 0.7);
    ASSERTEQUALT(2.338, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<float>(0.699, 1.643), 0.001);
    E = phf.Ekt(E, 0.8);
    ASSERTEQUALT(2.741, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<float>(0.904, 1.298), 0.001);
    E = phf.Ekt(E, 0.9);
    ASSERTEQUALT(3.015, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<float>(0.990, 1.038), 0.001);
    E = phf.Ekt(E, 1.0);
    ASSERTEQUALT(3.055, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<float>(1.000, 1.000), 0.001);

    cout << "TEST	:   test_PH5Feed() OK " << endl;
}

void test_SamplePHCurve() {
    cout << "TEST	: test_SamplePHCurve()" << endl;

    // Create a PHCurve interpolator from the PH curve coefficients
    // for an arc connecting the three points {(-1,1),(0,2),(1,1)}
    PH5Curve<float> ph(ph_arc());

    float vMax = 100;		// maximum velocity (mm/s)
    float tvMax = 0.01;		// time to achieve maximum velocity (seconds)
    float vIn = 0;			// initial velocity (mm/s)
    float vCruise = vMax;	// cruising velocity (mm/s)
    float vOut = 0;			// final velocity (mm/s)

    // Create a quintic feedrate for traversing the above curve smoothly
    PHFeed<float> phf(ph, vMax, tvMax, vIn, vCruise, vOut);

    int N = 100;			// number of points to interpolate
    float E = 0;			// interpolation normalized parametric state [0,1]

    // Generate a set of points using PH feed rate along PH curve
    cout << "X,Y" << endl;
    for (float tPoint = 0; tPoint <= N; tPoint++) {
        E = phf.Ekt(E, tPoint / N);
        Complex<float> point = ph.r(E);
        cout << point.Re() << "," << point.Im() << endl;
    }

    cout << "TEST	:   test_SamplePHCurve() OK " << endl;
}

int main(int argc, char *argv[]) {
    LOGINFO3("INFO	: ph5 test v%d.%d.%d", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH);
    firelog_level(FIRELOG_TRACE);

    test_Complex();
    test_Bernstein();
    test_PH5Curve();
    test_PHFeed();
    test_SamplePHCurve();

    cout << "TEST	: END OF TEST main()" << endl;
}
