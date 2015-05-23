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
template <class T>
class CLASS_DECLSPEC tinyvector {
    private:
        T elt[VECTOR_SIZE];
        int16_t length;

    public:
        tinyvector() : length(0) {
        }
        bool push_back(const T& val) {
            if (length >= VECTOR_SIZE) {
                return false;
            }
            elt[length++] = val;
            return true;
        }
        T& operator[](size_t index) {
            return elt[index];
        }
        const T& operator[](size_t index) const {
            return elt[index];
        }
        size_t size() {
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
        string stringify(int nPlaces = 0);
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
        int N;

    protected:
        Complex<T> calc_wij(int i, int j);
        Complex<T> pik(int i, int k);
        Complex<T> rit(int i, T p);
        T sit(int i, T p);
        T sik(int i, int j);
        T sigmaij(int i, int j);
        Complex<T> ritprime(int i, T p);

    public:
        PH5Curve(PHVECTOR<Complex<T> > phz, PHVECTOR<Complex<T> > phq);
        Complex<T> r(T p);
        T s(T p);
        T sigma(T p);
        Complex<T> rprime(T p);
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
        int iterations;		// maximum Newton-Raphson iterations
        T tauCruise;
        T tauDecel;
        T Faccel[7];
        T Fcruise[7];
        T Fdecel[7];
        inline T Vaccel(int k) {
            return k < 3 ? vIn : vCruise;
        };
        inline T Vdecel(int k) {
            return k < 3 ? vCruise : vOut;
        };
        T Fk(T vIn, T vOut, int k);
        inline T Vk(T vIn, T vOut, int k) {
            return k < 3 ? vIn : vOut;
        }

    public:
        PHFeed(PH5Curve<T> &ph5, T vMax = 200, T tvMax = 0.1, T vIn = 0, T vCruise = 200, T vOut = 0);
    public:
        T Ft(T tau);
    public:
        T Ekt(T Ekprev, T tau);

    public:
        inline T sigma(T E) {
            return ph.sigma(E);
        }
        inline T s(T E) {
            return ph.s(E);
        }
        inline Complex<T> r(T E) {
            return ph.r(E);
        }
        inline T get_tS() {
            return tS;
        };
        inline T get_tAccel() {
            return tAccel;
        };
        inline T get_tCruise() {
            return tCruise;
        };
        inline T get_tDecel() {
            return tDecel;
        };
        inline T get_sAccel() {
            return sAccel;
        };
        inline T get_sCruise() {
            return sCruise;
        };
        inline T get_sDecel() {
            return sDecel;
        };
        inline T get_S() {
            return S;
        };
        inline T get_vIn() {
            return vIn;
        };
        inline T get_vCruise() {
            return vCruise;
        };
        inline T get_vOut() {
            return vOut;
        };
        inline T get_vMax() {
            return vMax;
        };
        inline T get_tvMax() {
            return tvMax;
        };

};

extern int choose5[6];
extern int choose6[7];

template <class T>
T Bernstein5(int k, T p) {
    T result = choose5[k];
    T p1 = 1 - p;
    for (int i = 0; i < 5 - k; i++) {
        result *= p1;
    }
    for (int i = 0; i < k; i++) {
        result *= p;
    }
    return result;
}

template <class T>
T Bernstein6(int k, T p) {
    T result = choose6[k];
    T p1 = 1 - p;
    for (int i = 0; i < 6 - k; i++) {
        result *= p1;
    }
    for (int i = 0; i < k; i++) {
        result *= p;
    }
    return result;
}

} // namespace ph5

#endif
