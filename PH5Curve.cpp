#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "FireLog.h"
#include "FireUtils.hpp"
#include "version.h"
#include "ph5.hpp"

using namespace std;
using namespace ph5;

namespace ph5 {
int choose5[6] = { 1, 5, 10, 10, 5, 1 };
int choose6[7] = { 1, 6, 15, 20, 15, 6, 1 };
}

template<class T>
PH5Curve<T>::PH5Curve(vector<Complex<T> > phz, vector<Complex<T> > phq) {
	ASSERTEQUAL(phz.size(), phq.size());
	ASSERT(phz.size() > 2);
	phz[0].assertEqualT(Complex<T>());
	
	this->z = phz;
	this->q = phq;
	this->N = phz.size() - 1;

	wi0.push_back(Complex<T>());
	wi1.push_back(Complex<T>());
	wi2.push_back(Complex<T>());
	for (int i=1; i<=N; i++) {
		wi0.push_back(calc_wij(i,0));
		wi1.push_back(calc_wij(i,1));
		wi2.push_back(calc_wij(i,2));
	}
	pi0.push_back(Complex<T>());
	pi1.push_back(Complex<T>());
	pi2.push_back(Complex<T>());
	pi3.push_back(Complex<T>());
	pi4.push_back(Complex<T>());
	pi5.push_back(Complex<T>());
	for (int i=1; i<=N; i++) {
		pi0.push_back(pik(i,0));
		pi1.push_back(pik(i,1));
		pi2.push_back(pik(i,2));
		pi3.push_back(pik(i,3));
		pi4.push_back(pik(i,4));
		pi5.push_back(pik(i,5));
	}
	sigmai0.push_back(0);
	sigmai1.push_back(0);
	sigmai2.push_back(0);
	sigmai3.push_back(0);
	sigmai4.push_back(0);
	for (int i=1; i<=N; i++) {
		sigmai0.push_back(sigmaij(i,0));
		sigmai1.push_back(sigmaij(i,1));
		sigmai2.push_back(sigmaij(i,2));
		sigmai3.push_back(sigmaij(i,3));
		sigmai4.push_back(sigmaij(i,4));
	}
}

template<class T>
T PH5Curve<T>::sigma(T p) {
	return rprime(p).modulus();
}

template<class T>
Complex<T> PH5Curve<T>::rprime(T p) {
	ASSERT(0 <= p);
	ASSERT(p <= 1);
	T PN = p * N;
	int iPN = p == 1 ? N : (((int) PN)+1);
	return ritprime(iPN, PN-iPN+1)*N;
}

template<class T>
Complex<T> PH5Curve<T>::ritprime(int i, T p) {
	Complex<T> sum;
	T p1 = 1 - p;
	if (i == 1) {
		sum.add((z[1]*3-z[2])*p1*p1/2.0);
		sum.add(z[1]*2*p1*p);
		sum.add((z[1]+z[2])*p*p/2.0);
	} else if (i == N) {
		sum.add((z[N]+z[N-1])*p1*p1/2.0);
		sum.add(z[N]*2*p1*p);
		sum.add((z[N]*3-z[N-1])*p*p/2.0);
	} else {
		sum.add((z[i-1]+z[i])*p1*p1/2.0);
		sum.add(z[i]*2*p1*p);
		sum.add((z[i]+z[i+1])*p*p/2.0);
	}
	return sum * sum;
}

template<class T>
T PH5Curve<T>::s(T p) {
	ASSERT(0 <= p);
	ASSERT(p <= 1);
	T PN = p * N;
	int iPN = p == 1 ? N : (((int) PN) + 1);
	T sum = 0;
	for (int iSeg=1; iSeg<iPN; iSeg++) {
		sum += sit(iSeg, 1);
	}
	sum += sit(iPN, PN-iPN+1);
	return sum;
}

template<class T>
T PH5Curve<T>::sit(int i, T p) {
	T sum = 0;
	for (int k=0; k<=5; k++) {
		sum += sik(i,k) * Bernstein5(k, p);
	}
}

template<class T>
T PH5Curve<T>::sik(int i, int k) {
    T sum = 0;
    for (int j=0; j<=k-1; j++) {
        switch (j) {
        case 0:
            sum += sigmai0[i];
            break;
        case 1:
            sum += sigmai1[i];
            break;
        case 2:
            sum += sigmai2[i];
            break;
        case 3:
            sum += sigmai3[i];
            break;
        case 4:
            sum += sigmai4[i];
            break;
		default:
			ASSERTFAIL("si?");
			break;
        }
    }
	return sum/5;
}

template<class T>
T PH5Curve<T>::sigmaij(int i, int j) {
    T u0 = wi0[i].Re();
    T v0 = wi0[i].Im();
    T u1 = wi1[i].Re();
    T v1 = wi1[i].Im();
    T u2 = wi2[i].Re();
    T v2 = wi2[i].Im();
    switch (j) {
    case 0:
        return u0 * u0 + v0 * v0;
    case 1:
        return u0 * u1 + v0 * v1;
    case 2:
        return (2 / 3.0) * (u1 * u1 + v1 * v1) + (1 / 3.0) * (u0 * u2 + v0 * v2);
    case 3:
        return u1 * u2 + v1 * v2;
    case 4:
        return u2 * u2 + v2 * v2;
    default:
        ASSERTFAIL("sigmai?");
    }
}

template<class T>
Complex<T> PH5Curve<T>::r(T p) {
    T PN = p * N;
    int i = p < 1 ? ((int)PN) + 1 : PN;
    return rit(i, PN - i + 1);
}

template<class T>
Complex<T> PH5Curve<T>::rit(int i, T e) {
    Complex<T> sum;
    T e1 = 1 - e;
    T ek[6];
    T e1k[6];
    ek[0] = e1k[5] = 1;
    T Eprod = 1;
    T E1prod = 1;
    for (int k = 1; k <= 5; k++) {
        ek[k] = Eprod = e * Eprod;
        e1k[5 - k] = E1prod = e1 * E1prod;
    }
    sum.add(pi0[i] * (choose5[0]*e1k[0]*ek[0]));
    sum.add(pi1[i] * (choose5[1]*e1k[1]*ek[1]));
    sum.add(pi2[i] * (choose5[2]*e1k[2]*ek[2]));
    sum.add(pi3[i] * (choose5[3]*e1k[3]*ek[3]));
    sum.add(pi4[i] * (choose5[4]*e1k[4]*ek[4]));
    sum.add(pi5[i] * (choose5[5]*e1k[5]*ek[5]));
    return sum;
}

template<class T>
Complex<T> PH5Curve<T>::calc_wij(int i, int j) {
    if (i == 1) {
        switch (j) {
        case 0:
            return (z[1] * 3 - z[2]) / 2;
        case 1:
            return z[1];
        case 2:
            return (z[1] + z[2]) / 2;
        default:
            ASSERTFAIL("w1?");
        }
    }
    if (i == N) {
        switch (j) {
        case 0:
            return (z[N - 1] + z[N]) / 2;
        case 1:
            return z[N];
        case 2:
            return (z[N] * 3 - z[N - 1]) / 2;
        default:
            ASSERTFAIL("wN?");
        }
    }
    switch (j) {
    case 0:
        return (z[i - 1] + z[i]) / 2;
    case 1:
        return z[i];
    case 2:
        return (z[i] + z[i + 1]) / 2;
    default:
        ASSERTFAIL("calc_wij j");
    }
}
template<class T>
Complex<T> PH5Curve<T>::pik(int i, int k) {
    ASSERT(i > 0);
    ASSERT(i <= N);
    switch (k) {
    case 0:
        return q[i - 1];
    case 1:
        return pik(i, 0) + wi0[i] * wi0[i] / 5;
    case 2:
        return pik(i, 1) + wi0[i] * wi1[i] / 5;
    case 3:
        return pik(i, 2) + wi1[i] * wi1[i] * 2 / 15 + wi0[i] * wi2[i] / 15;
    case 4:
        return pik(i, 3) + wi1[i] * wi2[i] / 5;
    case 5:
        return pik(i, 4) + wi2[i] * wi2[i] / 5;
    default:
        ASSERTFAIL("invalid k");
    }
}

template class PH5Curve<float>;
template class PH5Curve<double>;

