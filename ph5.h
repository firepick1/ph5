#ifndef PH5_H
#define PH5_H

#ifdef _MSC_VER
#include "winjunk.hpp"
#else
#define CLASS_DECLSPEC
#endif

#ifdef CMAKE
#include <stdint.h>
#include <cstdio>
#include <string>
#endif

#ifdef ARDUINO
#include "Arduino.h"
#define MEMORY_MODEL_TINY
#endif

#ifdef MEMORY_MODEL_TINY
// MEMORY_MODEL_TINY: We have no room for niceties
#define PHVECTOR tinyvector
#define VECTOR_SIZE 3
#define PH5TYPE float
#else
// MEMORY_MODEL_LARGE: do it right
#include <vector>
#include <cstdlib>
#include "FireLog.h"
#include "FireUtils.hpp"
#define PHVECTOR vector
#define PH5TYPE double
#endif

#ifndef ASSERT
#define ASSERT(x)
#define ASSERTFAIL(x)
#define ASSERTEQUALT(a,b,c)
#define ASSERTEQUAL(a,b)
#define ASSERTEQUALS(a,b)
#define firelog_level(x)
#define LOGERROR1(a,b)
#define LOGDEBUG1(a,b)
#define LOGINFO1(a,b)
#define LOGINFO2(a,b,c)
#define LOGINFO3(a,b,c,d)
#endif

#include <math.h>

using namespace std;

#if __amd64__ || __x86_64__ || _WIN64 || _M_X64
#define PH5_64_BIT
#define PH5_PLATFORM_BITS 64
#else
#define PH5_32_BIT
#define PH5_PLATFORM_BITS 32
#endif

namespace ph5 {

#ifdef MEMORY_MODEL_TINY
extern PH5TYPE OUTOFRANGE;

template <class T>
class CLASS_DECLSPEC tinyvector {
    private:
        T elt[VECTOR_SIZE];
        int16_t length;

    public:
        tinyvector() : length(0) {
        }
        inline bool push_back(const T& val) {
            if (length >= VECTOR_SIZE) {
#ifdef CMAKE
				throw "HELP";
#endif
                return false;
            }
            elt[length++] = val;
            return true;
        }
        inline T& operator[](int16_t index) {
            return elt[index];
        } 
		inline const T& operator[](int16_t index) const {
			if (index < 0 || VECTOR_SIZE <= index) {
				return OUTOFRANGE;
			}
            return elt[index];
        }
        inline size_t size() {
            return length;
        }
};
#endif

template <class T>
class CLASS_DECLSPEC Complex {
    private:
        T re;
        T im;

    public:
        Complex(T re = 0, T im = 0);
    public:
        inline T Re() {
            return this->re;
        }
    public:
        inline T Im() {
            return this->im;
        }
    public:
        T modulus() {
            return ::sqrt(re * re + im * im);
        }
    public:
        inline Complex<T> operator+(Complex<T> that) {
            return Complex(re + that.re, im + that.im);
        }
    public:
        inline Complex<T> operator+(T k) {
            return Complex(re + k, im);
        }
    public:
        inline Complex<T> operator-(Complex<T> that) {
            return Complex(re - that.re, im - that.im);
        }
    public:
        inline Complex<T> operator-(T k) {
            return Complex(re - k, im);
        }
    public:
        inline bool operator==(Complex<T> that) {
            return re == that.re && im == that.im;
        }
    public:
        inline Complex<T> operator*(Complex<T> that) {
            return Complex(re * that.re - im * that.im, re * that.im + im * that.re);
        }
    public:
        inline Complex<T> operator*(T that) {
            return Complex(re * that, im * that);
        }
    public:
        inline Complex<T> conj() {
            return Complex(re, -im);
        }
    public:
        inline Complex<T> recip() {
            T denom = re * re + im * im;
            return Complex(re / denom, -im / denom);
        }
    public:
        inline Complex<T> operator/ (Complex<T> c2) {
            T denom = c2.re * c2.re + c2.im * c2.im;
            return Complex(
                       (re * c2.re + im * c2.im) / denom,
                       (im * c2.re - re * c2.im) / denom
                   );
        }
    public:
        inline Complex<T> operator/ (T k) {
            return Complex(re / k, im / k);
        }
    public:
        inline void add(Complex<T> that) {
            re += that.re;
            im += that.im;
        }
    public:
        inline Complex<T> sqrt() {
            T m = modulus();
            T p = ::sqrt((m + re) / 2);
            T q = ::sqrt((m - re) / 2);
            if (im >= 0) {
                return Complex(p, q);
            } else {
                return Complex(p, -q);
            }
        }
#ifdef CMAKE
    public:
        string stringify(int16_t nPlaces = 0);
#endif

    public:
        bool assertEqualT(Complex<T> that, double tolerance = 0.0000001);
};

namespace ph5 {
template<class T>
inline Complex<T> operator*(T k, Complex<T> c) {
    return Complex<T>(k * c.re, k * c.im);
}
}

template<class T>
class CLASS_DECLSPEC PH5Curve {
    private:
        PHVECTOR<Complex<T> > z;
        PHVECTOR<Complex<T> > q;
        PHVECTOR<Complex<T> > wi0;
        PHVECTOR<Complex<T> > wi1;
        PHVECTOR<Complex<T> > wi2;
        PHVECTOR<Complex<T> > pi0;
        PHVECTOR<Complex<T> > pi1;
        PHVECTOR<Complex<T> > pi2;
        PHVECTOR<Complex<T> > pi3;
        PHVECTOR<Complex<T> > pi4;
        PHVECTOR<Complex<T> > pi5;
        PHVECTOR<T> sigmai0;
        PHVECTOR<T> sigmai1;
        PHVECTOR<T> sigmai2;
        PHVECTOR<T> sigmai3;
        PHVECTOR<T> sigmai4;
        int16_t N;

    protected:
        Complex<T> calc_wij(int16_t i, int16_t j);
        Complex<T> pik(int16_t i, int16_t k);
        Complex<T> rit(int16_t i, T p);
        T sit(int16_t i, T p);
        T sik(int16_t i, int16_t j);
        T sigmaij(int16_t i, int16_t j);
        Complex<T> ritprime(int16_t i, T p);

    public:
        PH5Curve(PHVECTOR<Complex<T> > phz, PHVECTOR<Complex<T> > phq);
        Complex<T> r(T p);
        T s(T p);
        T sigma(T p);
        Complex<T> rprime(T p);
};

enum FeedUseCase {
	FEEDUSE_A,	// vIn==vCruise==vOut > 0
	FEEDUSE_B1,	// vOut==vCruise sAccel<sMax
	FEEDUSE_B2,	// vOut==vCruise sAccel:sMax
	FEEDUSE_C1,	// vIn==vCruise sAccel<sMax
	FEEDUSE_C2,	// vIn==vCruise sAccel:sMax
	FEEDUSE_D1, // vIn:0, vOut:0, vCruise<vMax, tCruise:0
	FEEDUSE_D2, // vIn:0, vOut:0, vCruise:vMax, tCruise>=0
};

template<class T>
class CLASS_DECLSPEC PHFeed {
    protected:
        PH5Curve<T> &ph;	// REQUIRED: curve to traverse
        T vMax;				// OPTION: maximum velocity
        T tvMax;			// OPTION: time from rest to maximum velocity
        T vIn;				// OPTION: PH curve entry velocity
        T vCruise;			// OPTION: PH curve constant cruise velocity
        T vOut;				// OPTION: PH curve exit velocity
        T S;				// curve arc length
        T sAccel;			// acceleration arc length distance
        T sCruise;			// cruise arc length distance
        T sDecel;			// deceleration arc length distance
        T tAccel;			// acceleration time
        T tCruise;			// cruise time
        T tDecel;			// deceleration time
        T tS;				// total traversal time
        T epsilon;			// Newton-Raphson convergence limit
        int16_t iterations;	// maximum Newton-Raphson iterations
        T tauCruise;
        T tauDecel;
        T Faccel[7];
        T Fcruise[7];
        T Fdecel[7];
		int8_t uc;			// FeedUseCase
        inline T Vaccel(int16_t k) {
            return k < 3 ? vIn : vCruise;
        };
        inline T Vdecel(int16_t k) {
            return k < 3 ? vCruise : vOut;
        };
        T Fk(T vIn, T vOut, int16_t k);
        inline T Vk(T vIn, T vOut, int16_t k) {
            return k < 3 ? vIn : vOut;
        }

    public:
        PHFeed(PH5Curve<T> &ph5, 
			T vMax = 200,	// maximum velocity 
			T tvMax = 0.1,	// time to reach maximum velocity 
			T vIn = 0,		// entry velocity 
			T vOut = 0,		// exit velocity
			T vCruise = 0 	// target cruise velocity (default is vMax)
			);
    public:
        T Ft(T tau);
    public:
        T Ekt(T Ekprev, T tau);

    public:
		FeedUseCase getFeedUseCase() {
			return (FeedUseCase) uc;
		}
        inline T sigma(T E) { // parametric arc length velocity
            return ph.sigma(E);
        }
        inline T s(T E) { // parametric distance
            return ph.s(E);
        }
        inline Complex<T> r(T E) { 
            return ph.r(E);
        }
        inline T get_tS() { // total traversal time
            return tS;
        };
        inline T get_tAccel() { // acceleration time
            return tAccel;
        };
        inline T get_tCruise() { // cruising time
            return tCruise;
        };
        inline T get_tDecel() { // deceleration time
            return tDecel;
        };
        inline T get_sAccel() { // acceleration distance
            return sAccel;
        };
        inline T get_sCruise() { // cruising distance
            return sCruise;
        };
        inline T get_sDecel() { // deceleration distance
            return sDecel;
        };
        inline T get_S() { // total path length
            return S;
        };
        inline T get_vIn() { // entry velocity
            return vIn;
        };
        inline T get_vCruise() { // cruise velocity
            return vCruise;
        };
        inline T get_vOut() { // exit velocity
            return vOut;
        };
        inline T get_vMax() { // maximum allowed velocity
            return vMax;
        };
        inline T get_tvMax() { // time to reach maximum allows velocity
            return tvMax;
        };

};

extern int16_t choose5[6];
extern int16_t choose6[7];

template <class T>
T Bernstein5(int16_t k, T p) {
    T result = choose5[k];
    T p1 = 1 - p;
    for (int16_t i = 0; i < 5 - k; i++) {
        result *= p1;
    }
    for (int16_t i = 0; i < k; i++) {
        result *= p;
    }
    return result;
}

template <class T>
T Bernstein6(int16_t k, T p) {
    T result = choose6[k];
    T p1 = 1 - p;
    for (int16_t i = 0; i < 6 - k; i++) {
        result *= p1;
    }
    for (int16_t i = 0; i < k; i++) {
        result *= p;
    }
    return result;
}

} // namespace ph5

#endif
