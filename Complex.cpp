#ifdef CMAKE
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#endif

#include <math.h>
#include "version.h"
#include "ph5.h"

using namespace std;
using namespace ph5;

template<class T>
Complex<T>::Complex(T re, T im) 
	: re(re), im(im) {
}

template<class T>
bool Complex<T>::assertEqualT(Complex<T> that, double tolerance) {
	ASSERTEQUALT(that.re, re, tolerance);
	ASSERTEQUALT(that.im, im, tolerance);
}

#ifdef CMAKE
template<class T>
string Complex<T>::stringify(int16_t nPlaces) {
	string s;
	char fmt[100];
	snprintf(fmt, sizeof(fmt), "%%.%dlf", nPlaces);
	char sre[100];
	snprintf(sre, sizeof(sre), fmt, (double) re);
	char sim[100];
	snprintf(sim, sizeof(sim), fmt, (double) im);
	if (im) {
		if (im < 0) {
			if (re != 0) {
				s.append(sre);
			}
			if (im == -1) {
				s.append("-");
			} else {
				s.append(sim);
			}
		} else {
			if (re != 0) {
				s.append(sre);
				s.append("+");
			}
			if (im != 1) {
				s.append(sim);
			}
		}
		s.append("i");
	} else {
		s.append(sre);
	}
	return s;
}
#endif

template class Complex<PH5TYPE>;
