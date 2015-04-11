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

#ifdef JS
})
#endif

void test_Complex() {
	Complex<float> c1(1,2);
	ASSERTEQUAL(c1.Re(), 1);
	ASSERTEQUAL(c1.Im(), 2);

	Complex<float> c5(3,4);
	ASSERTEQUAL(5, c5.modulus());

	Complex<float> c15 = c1+c5;
	ASSERT(c15 == Complex<float>(4,6));
	ASSERT(c1+1 == Complex<float>(2,2));

	Complex<float> c13(1,3);
	Complex<float> c24(2,4);
	Complex<float> c1324 = c13-c24;
	ASSERT((c13-c24) == Complex<float>(-1,-1));
	ASSERT((c24-c13) == Complex<float>(1,1));
	ASSERT((c13-1) == Complex<float>(0,3));

	ASSERT((c13*c24) == Complex<float>(1*2-3*4,1*4+3*2));
	ASSERT((c24*c13) == Complex<float>(1*2-3*4,1*4+3*2));
	ASSERT((c13*2) == Complex<float>(2,6));

	ASSERT(c13.conj() == Complex<float>(1,-3));

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
	ASSERT(sum1 == Complex<float>(1,3));
	sum1.add(c24);
	ASSERT(sum1 == Complex<float>(3,7));

	ASSERT(Complex<float>(25).sqrt() == Complex<float>(5));
	ASSERT(Complex<float>(-1).sqrt() == Complex<float>(0,1));
	ASSERT(Complex<float>(3,4).sqrt() == Complex<float>(2,1));

	ASSERTEQUALS("0", Complex<float>().stringify().c_str());
	ASSERTEQUALS("1", Complex<float>(1).stringify().c_str());
	ASSERTEQUALS("i", Complex<float>(0,1).stringify().c_str());
	ASSERTEQUALS("-i", Complex<float>(0,-1).stringify().c_str());
	ASSERTEQUALS("2i", Complex<float>(0,2).stringify().c_str());
	ASSERTEQUALS("-2i", Complex<float>(0,-2).stringify().c_str());
	ASSERTEQUALS("1+i", Complex<float>(1,1).stringify().c_str());
	ASSERTEQUALS("1-i", Complex<float>(1,-1).stringify().c_str());
	ASSERTEQUALS("1+2i", Complex<float>(1,2).stringify().c_str());
	ASSERTEQUALS("1-2i", Complex<float>(1,-2).stringify().c_str());
	ASSERTEQUALS("1.2-5.7i", Complex<float>(1.2345,-5.6789).stringify(1).c_str());
	ASSERTEQUALS("1.0+1.0i", Complex<float>(0.99,0.99).stringify(1).c_str());
	ASSERTEQUALS("1.0-1.0i", Complex<float>(0.99,-0.99).stringify(1).c_str());
	ASSERTEQUALS("0.99-0.99i", Complex<float>(0.99,-0.99).stringify(2).c_str());
}

int main(int argc, char *argv[])
{
    LOGINFO3("ph5 test v%d.%d.%d", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH);
    firelog_level(FIRELOG_TRACE);

    cout << "test_Complex()" << endl;
    test_Complex();

    cout << "END OF TEST main()" << endl;
}
