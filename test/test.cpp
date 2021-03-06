#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "FireLog.h"
#include "FireUtils.hpp"
#include "ph5.h"
#include "version.h"

using namespace std;
using namespace ph5;

void test_Complex() {
    cout << "TEST	: test_Complex()" << endl;

    Complex<PH5TYPE> c1(1, 2);
    ASSERTEQUAL(c1.Re(), 1);
    ASSERTEQUAL(c1.Im(), 2);

    Complex<PH5TYPE> c5(3, 4);
    ASSERTEQUAL(5, c5.modulus());

    Complex<PH5TYPE> c15 = c1 + c5;
    ASSERT(c15 == Complex<PH5TYPE>(4, 6));
    ASSERT(c1 + 1 == Complex<PH5TYPE>(2, 2));

    Complex<PH5TYPE> c13(1, 3);
    Complex<PH5TYPE> c24(2, 4);
    Complex<PH5TYPE> c1324 = c13 - c24;
    ASSERT((c13 - c24) == Complex<PH5TYPE>(-1, -1));
    ASSERT((c24 - c13) == Complex<PH5TYPE>(1, 1));
    ASSERT((c13 - 1) == Complex<PH5TYPE>(0, 3));

    ASSERT((c13 * c24) == Complex<PH5TYPE>(1 * 2 - 3 * 4, 1 * 4 + 3 * 2));
    ASSERT((c24 * c13) == Complex<PH5TYPE>(1 * 2 - 3 * 4, 1 * 4 + 3 * 2));
    ASSERT((c13 * 2) == Complex<PH5TYPE>(2, 6));

    ASSERT(c13.conj() == Complex<PH5TYPE>(1, -3));

    Complex<PH5TYPE> cr = c13.recip();
    PH5TYPE epsilon = 0.000001;
    ASSERTEQUALT(cr.Re(), 0.1, epsilon);
    ASSERTEQUALT(cr.Im(), -0.3, epsilon);
    Complex<PH5TYPE> crr = cr.recip();
    ASSERTEQUALT(crr.Re(), 1, epsilon);
    ASSERTEQUALT(crr.Im(), 3, epsilon);

    Complex<PH5TYPE> c13div24 = c13 / c24;
    ASSERTEQUALT(c13div24.Re(), 0.7, epsilon);
    ASSERTEQUALT(c13div24.Im(), 0.1, epsilon);

    Complex<PH5TYPE> sum1;
    sum1.add(c13);
    ASSERT(sum1 == Complex<PH5TYPE>(1, 3));
    sum1.add(c24);
    ASSERT(sum1 == Complex<PH5TYPE>(3, 7));

    ASSERT(Complex<PH5TYPE>(25).sqrt() == Complex<PH5TYPE>(5));
    ASSERT(Complex<PH5TYPE>(-1).sqrt() == Complex<PH5TYPE>(0, 1));
    ASSERT(Complex<PH5TYPE>(3, 4).sqrt() == Complex<PH5TYPE>(2, 1));

    ASSERTEQUALS("0", Complex<PH5TYPE>().stringify().c_str());
    ASSERTEQUALS("1", Complex<PH5TYPE>(1).stringify().c_str());
    ASSERTEQUALS("i", Complex<PH5TYPE>(0, 1).stringify().c_str());
    ASSERTEQUALS("-i", Complex<PH5TYPE>(0, -1).stringify().c_str());
    ASSERTEQUALS("2i", Complex<PH5TYPE>(0, 2).stringify().c_str());
    ASSERTEQUALS("-2i", Complex<PH5TYPE>(0, -2).stringify().c_str());
    ASSERTEQUALS("1+i", Complex<PH5TYPE>(1, 1).stringify().c_str());
    ASSERTEQUALS("1-i", Complex<PH5TYPE>(1, -1).stringify().c_str());
    ASSERTEQUALS("1+2i", Complex<PH5TYPE>(1, 2).stringify().c_str());
    ASSERTEQUALS("1-2i", Complex<PH5TYPE>(1, -2).stringify().c_str());
    ASSERTEQUALS("1.2-5.7i", Complex<PH5TYPE>(1.2345, -5.6789).stringify(1).c_str());
    ASSERTEQUALS("1.0+1.0i", Complex<PH5TYPE>(0.99, 0.99).stringify(1).c_str());
    ASSERTEQUALS("1.0-1.0i", Complex<PH5TYPE>(0.99, -0.99).stringify(1).c_str());
    ASSERTEQUALS("0.99-0.99i", Complex<PH5TYPE>(0.99, -0.99).stringify(2).c_str());

    c13.assertEqualT(Complex<PH5TYPE>(1.0009, 2.9991), 0.001);
    cout << "TEST	:   test_Complex() OK" << endl;
}

void test_Bernstein() {
    PH5TYPE  epsilon = 0.00001;
    cout << "TEST	: test_Bernstein()" << endl;
    ASSERTEQUALT(0, Bernstein5<PH5TYPE>(5, 0), epsilon);
    ASSERTEQUALT(0.03125, Bernstein5<PH5TYPE>(5, 0.5), epsilon);
    ASSERTEQUALT(0.15625, Bernstein5<PH5TYPE>(1, 0.5), epsilon);
    cout << "TEST	:   test_Bernstein() OK" << endl;
}

PH5Curve<PH5TYPE> ph_arc() {
    PHVECTOR<Complex<PH5TYPE> > q;
    PHVECTOR<Complex<PH5TYPE> > z;
    q.push_back(Complex<PH5TYPE>(-1, 1));
    q.push_back(Complex<PH5TYPE>(0, 2));
    q.push_back(Complex<PH5TYPE>(1, 1));
    z.push_back(Complex<PH5TYPE>());
    z.push_back(Complex<PH5TYPE>(1.124171968973597, 0.444771808762066));
    z.push_back(Complex<PH5TYPE>(1.124171968973597, -0.444771808762066));
    return PH5Curve<PH5TYPE>(z, q);
}

void test_PH5Curve() {
    cout << "TEST	: test_PH5Curve()" << endl;
    PH5Curve<PH5TYPE> ph(ph_arc());

    long msStart = millis();
    PH5TYPE epsilon = 0.00001;
    int16_t ITER = 10000;
    for (int16_t i = 0; i < ITER; i++) {
        Complex<PH5TYPE> c;
        c = ph.r(0);
        c.assertEqualT(Complex<PH5TYPE>(-1, 1), epsilon);
        c = ph.r(0.25);
        c.assertEqualT(Complex<PH5TYPE>(-0.598911, 1.75), epsilon);
        c = ph.r(0.5);
        c.assertEqualT(Complex<PH5TYPE>(0, 2), epsilon);
        c = ph.r(0.75);
        c.assertEqualT(Complex<PH5TYPE>(0.598911, 1.75), epsilon);
        c = ph.r(1.0);
        c.assertEqualT(Complex<PH5TYPE>(1, 1), epsilon);
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
    ph.rprime(0.00).assertEqualT(Complex<PH5TYPE>(0.945, 4.000), epsilon);
    ph.rprime(0.25).assertEqualT(Complex<PH5TYPE>(2.132, 2.000), epsilon);
    ph.rprime(0.50).assertEqualT(Complex<PH5TYPE>(2.528, 0.000), epsilon);
    ph.rprime(0.75).assertEqualT(Complex<PH5TYPE>(2.132, -2.000), epsilon);
    ph.rprime(1.00).assertEqualT(Complex<PH5TYPE>(0.945, -4.000), epsilon);

    cout << "TEST	:   test_PH5Curve() OK " << endl;
}

void test_PHFeed() {
    cout << "TEST	: test_PH5Feed()" << endl;
    PH5Curve<PH5TYPE> ph(ph_arc());
    PHFeed<PH5TYPE> phf(ph, 100, 0.01);
    PH5TYPE E = 0;
    PH5TYPE epsilon = 0.001;
    ASSERTEQUALT(0, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<PH5TYPE>(-1.000, 1.000), epsilon);
    E = phf.Ekt(E, 0.1);
    ASSERTEQUALT(0.039, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<PH5TYPE>(-0.990, 1.038), 0.001);
    E = phf.Ekt(E, 0.2);
    ASSERTEQUALT(0.313, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<PH5TYPE>(-0.904, 1.298), 0.001);
    E = phf.Ekt(E, 0.3);
    ASSERTEQUALT(0.716, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<PH5TYPE>(-0.699, 1.643), 0.001);
    E = phf.Ekt(E, 0.4);
    ASSERTEQUALT(1.122, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<PH5TYPE>(-0.389, 1.901), 0.001);
    E = phf.Ekt(E, 0.5);
    ASSERTEQUALT(1.527, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<PH5TYPE>(-0.000, 2.000), 0.001);
    E = phf.Ekt(E, 0.6);
    ASSERTEQUALT(1.933, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<PH5TYPE>(0.389, 1.901), 0.001);
    E = phf.Ekt(E, 0.7);
    ASSERTEQUALT(2.338, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<PH5TYPE>(0.699, 1.643), 0.001);
    E = phf.Ekt(E, 0.8);
    ASSERTEQUALT(2.741, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<PH5TYPE>(0.904, 1.298), 0.001);
    E = phf.Ekt(E, 0.9);
    ASSERTEQUALT(3.015, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<PH5TYPE>(0.990, 1.038), 0.001);
    E = phf.Ekt(E, 1.0);
    ASSERTEQUALT(3.055, ph.s(E), epsilon);
    ph.r(E).assertEqualT(Complex<PH5TYPE>(1.000, 1.000), 0.001);

    cout << "TEST	:   test_PH5Feed() OK " << endl;
}

void test_SamplePHCurve() {
    cout << "TEST	: test_SamplePHCurve()" << endl;

    // Create a PHCurve interpolator from the PH curve coefficients
    // for an arc connecting the three points {(-1,1),(0,2),(1,1)}
    PH5Curve<PH5TYPE> ph(ph_arc());

    PH5TYPE vMax = 100; // maximum velocity (mm/s)
    PH5TYPE tvMax = 0.01; // time to achieve maximum velocity (seconds)
    PH5TYPE vIn = 0; // initial velocity (mm/s)
    PH5TYPE vCruise = vMax; // cruising velocity (mm/s)
    PH5TYPE vOut = 0; // final velocity (mm/s)

    // Create a quintic feedrate for traversing the above curve smoothly
    PHFeed<PH5TYPE> phf(ph, vMax, tvMax, vIn, vOut, vCruise);

    int16_t N = 100; // number of points to interpolate
    PH5TYPE E = 0; // interpolation normalized parametric state [0,1]

    // Generate a set of points using PH feed rate along PH curve
    cout << "X,Y" << endl;
    for (PH5TYPE tPoint = 0; tPoint <= N; tPoint++) {
        E = phf.Ekt(E, tPoint / N);
        Complex<PH5TYPE> point = ph.r(E);
        cout << point.Re() << "," << point.Im() << endl;
    }

    cout << "TEST	:   test_SamplePHCurve() OK " << endl;
}

void test_size() {
	cout << "sizeof(Complex):" << sizeof(Complex<PH5TYPE>) << endl;
	cout << "sizeof(PH5Curve):" << sizeof(PH5Curve<PH5TYPE>) << endl;
	cout << "sizeof(PHFeed):" << sizeof(PHFeed<PH5TYPE>) << endl;
}

void test_performance() {
    PH5Curve<PH5TYPE> ph(ph_arc());
    PHFeed<PH5TYPE> phf(ph, 100, 0.01);
    int16_t N = 100; // number of points to interpolate
    PH5TYPE E = 0; // interpolation normalized parametric state [0,1]
	int iterations = 1000000;

	int32_t msStart;
	int32_t msElapsed;

	msStart = millis();
	PH5TYPE s;
	PH5TYPE Ekr = 0.35; // arbitrary
	for (int iter=0; iter<iterations; iter++) {
		s = ph.s(Ekr); // 269-330ms baseline;  268-330 improved
	}
	msElapsed = millis() - msStart;
	cout << "testPerformance ph.s():" << s << " " << msElapsed << "ms"<< endl;

	msStart = millis();
	PH5TYPE sigma;
	for (int iter=0; iter<iterations; iter++) {
		sigma = ph.sigma(Ekr); // 309-330ms baseline; 205-213 improved
	}
	msElapsed = millis() - msStart;
	cout << "testPerformance ph.sigma():" << sigma << " " << msElapsed << "ms"<< endl;

	msStart = millis();
	Complex<PH5TYPE> r;
	for (int iter=0; iter<iterations; iter++) {
		r = ph.r(Ekr); // 196-210ms baseline
	}
	msElapsed = millis() - msStart;
	cout << "testPerformance ph.r():" << msElapsed << "ms"<< endl;

	msStart = millis();
	PH5TYPE tau = 0.81;
	PH5TYPE Ftau;
	for (int iter=0; iter<iterations; iter++) {
		Ftau = phf.Ft(tau); // 185-205ms baseline
	}
	msElapsed = millis() - msStart;
	cout << "testPerformance phf.Ft():" << msElapsed << "ms"<< endl;

	msStart = millis();
	for (int iter=0; iter<1000; iter++) {
		for (PH5TYPE tPoint = 0; tPoint <= N; tPoint++) {
			E = phf.Ekt(E, tPoint / N);
			Complex<PH5TYPE> point = ph.r(E);
		}
	}
	msElapsed = millis() - msStart; // 231-240ms improved
	cout << "testPerformance phf.Ekt() + ph.r(E):" << msElapsed << endl;
}

int main(int argc, char *argv[]) {
    LOGINFO3("INFO	: ph5 test v%d.%d.%d", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH);
    firelog_level(FIRELOG_TRACE);

    test_Complex();
    test_Bernstein();
    test_PH5Curve();
    test_PHFeed();
    test_SamplePHCurve();
	test_size();
	test_performance();

    cout << "TEST	: END OF TEST main()" << endl;
}
