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

#ifdef JS
(function(firepick) {
    var DEGREE = 5;
    var b4 = new Bernstein(4);
    var b5 = new Bernstein(5);
    function PH5Curve(phz, phq, options) {
        var that = this;
        phq.length.should.equal(phz.length);
        options = options || {};
        that.logger = options.logger || new Logger(options);
        for (var i = 1; i < phz.length; i++) {
            that.logger.trace("phz[", i, "]:", phz[i]);
            phz[i].re.should.be.Number;
            phz[i].im.should.be.Number;
            phq[i].re.should.be.Number;
            phq[i].im.should.be.Number;
        }
        that.N = phz.length - 1;
        that.z = phz;
        that.q = phq;
        return that;
    };

    /////////////// PRIVATE ////////////////
    function powert(p, tk, t1k, K) {
        var p1 = 1 - p;
        tk.push(1);
        t1k.push(1);
        for (var k = 1; k <= K; k++) {
            tk.push(p * tk[k - 1]);
            t1k.splice(0, 0, p1 * t1k[0]);
        }
    };

    ///////////////// INSTANCE ///////////////
    PH5Curve.prototype.s = function(p) { // arc length
        var that = this;
        p.should.not.be.below(0);
        p.should.not.be.above(1);
        var PN = p * that.N;
        var iPN = Math.ceil(PN) || 1;
        var sum = 0;
        for (var iSeg = 1; iSeg < iPN; iSeg++) {
            sum += that.sit(iSeg, 1);
        }
        sum += that.sit(iPN, PN - iPN + 1);
        return sum;
    };
    PH5Curve.prototype.sit = function(i, p) { // arc length
        var that = this;
        var sum = 0;
        for (var k = 0; k <= DEGREE; k++) {
            var b5c = b5.coefficient(k, p);
            sum += that.sik(i, k) * b5c;
            that.logger.trace("sit k:", k, " sum:", sum, " b5c:", b5c, " p:", p);
        }
        return sum;
    };
    PH5Curve.prototype.sik = function(i, k) { // arc length
        var that = this;
        var sum = 0;
        for (var j = 0; j <= k - 1; j++) {
            sum += that.sigmaij(i, j);
        }
        return sum / DEGREE;
    };
    PH5Curve.prototype.sigmaij = function(i, j) {
        var that = this;
        var wi0 = that.wij(i, 0);
        var wi1 = that.wij(i, 1);
        var wi2 = that.wij(i, 2);
        var u0 = wi0.re;
        var v0 = wi0.im;
        var u1 = wi1.re;
        var v1 = wi1.im;
        var u2 = wi2.re;
        var v2 = wi2.im;
        switch (j) {
        case 0:
            return u0 * u0 + v0 * v0;
        case 1:
            return u0 * u1 + v0 * v1;
        case 2:
            return (2 / 3) * (u1 * u1 + v1 * v1) + (1 / 3) * (u0 * u2 + v0 * v2);
        case 3:
            return u1 * u2 + v1 * v2;
        case 4:
            return u2 * u2 + v2 * v2;
        default:
            should.fail("invalid j:" + j);
        }
    };
    PH5Curve.prototype.sigma = function(p) { // curve parametric speed
        var that = this;
        return that.rprime(p).modulus();
    };
    PH5Curve.prototype.rprime = function(p) { // hodograph
        var that = this;
        p.should.not.be.below(0);
        p.should.not.be.above(1);
        var PN = p * that.N;
        var i = Math.ceil(PN) || 1;
        return Complex.times(that.N, that.ritprime(i, PN - i + 1));
    };
    PH5Curve.prototype.ritprime = function(i, p) { // segment hodograph
        var that = this;
        var sum = new Complex();
        var p1 = 1 - p;
        var z = that.z;
        var N = that.N;
        if (i == = 1) {
            var z1 = z[1];
            var z2 = z[2];
            sum.add(Complex.times(1 / 2 * p1 * p1, Complex.times(3, z1).minus(z2)));
            sum.add(Complex.times(2 * p1 * p, z1));
            sum.add(Complex.times(1 / 2 * p * p, z1.plus(z2)));
        } else if (i == = N) {
            var zN = z[N];
            var zN1 = z[N - 1];
            sum.add(Complex.times(1 / 2 * p1 * p1, zN.plus(zN1)));
            sum.add(Complex.times(2 * p1 * p, zN));
            sum.add(Complex.times(1 / 2 * p * p, Complex.times(3, zN).minus(zN1)));
        } else {
            sum.add(Complex.times(1 / 2 * p1 * p1, z[i - 1].plus(z[i])));
            sum.add(Complex.times(2 * p1 * p, z[i]));
            sum.add(Complex.times(1 / 2 * p * p, z[i].plus(z[i + 1])));
        }
        return sum.times(sum);
    };
    PH5Curve.prototype.r = function(p) {
        var that = this;
        p.should.not.be.below(0);
        p.should.not.be.above(1);
        var PN = p * that.N;
        var i = Math.ceil(PN) || 1;
        return that.rit(i, PN - i + 1);
    };
    PH5Curve.prototype.rit = function(i, p) {
        var that = this;
        i.should.not.be.below(0);
        i.should.not.be.above(that.N);
        p.should.not.be.below(0);
        p.should.not.be.above(1);
        that.logger.trace("rit(", i, ",", p, ")");
        var sum = new Complex();
        var tk = [];
        var t1k = [];
        powert(p, tk, t1k, 5);
        for (var k = 0; k <= 5; k++) {
            var re = Util.choose(5, k) * t1k[k] * tk[k];
            var c = Complex.times(that.pik(i, k), re);
            sum.add(c);
            that.logger.trace("rit k:", k, " re:", re, " c:", c, " sum:", sum,
                              " pik:", that.pik(i, k), " choose:", Util.choose(5, k));
        }
        return sum;
    };
    PH5Curve.prototype.w1j = function(j) {
        var that = this;
        var z1 = that.z[1];
        var z2 = that.z[2];
        switch (j) {
        case 0:
            return Complex.times(1 / 2, Complex.times(3, z1).minus(z2));
        case 1:
            return z1;
        case 2:
            return Complex.times(1 / 2, z1.plus(z2));
        default:
            should.fail("w1j j:" + j);
        }
    };
    PH5Curve.prototype.wNj = function(j) {
        var that = this;
        var zN = that.z[that.N];
        var zN1 = that.z[that.N - 1];
        switch (j) {
        case 0:
            return Complex.times(1 / 2, zN1.plus(zN));
        case 1:
            return zN;
        case 2:
            return Complex.times(1 / 2, Complex.times(3, zN).minus(zN1));
        default:
            should.fail("wNj j:" + j);
        }
    };
    PH5Curve.prototype.wij = function(i, j) {
        var that = this;
        if (i == = 1) {
            return that.w1j(j);
        }
        if (i == = that.N) {
            return that.wNj(j);
        }
        var zi = that.z[i];
        i.should.not.be.below(1);
        i.should.not.be.above(that.N);
        zi.should.instanceOf(Complex);
        that.z[i - 1].should.instanceOf(Complex);
        switch (j) {
        case 0:
            return Complex.times(1 / 2, that.z[i - 1].plus(zi));
        case 1:
            return zi;
        case 2:
            return Complex.times(1 / 2, zi.plus(that.z[i + 1]));
        default:
            should.fail("wij j:" + j);
        }
    };
    PH5Curve.prototype.pik = function(i, k) {
        var that = this;
        i.should.be.above(0);
        i.should.not.be.above(that.N);

        switch (k) {
        case 0:
            return that.q[i - 1];
        case 1:
            return that.pik(i, 0)
                   .plus(Complex.times(1 / 5, that.wij(i, 0).times(that.wij(i, 0))));
        case 2:
            return that.pik(i, 1)
                   .plus(Complex.times(1 / 5, that.wij(i, 0).times(that.wij(i, 1))));
        case 3:
            return that.pik(i, 2)
                   .plus(Complex.times(2 / 15, that.wij(i, 1).times(that.wij(i, 1))))
                   .plus(Complex.times(1 / 15, that.wij(i, 0).times(that.wij(i, 2))));
        case 4:
            return that.pik(i, 3)
                   .plus(Complex.times(1 / 5, that.wij(i, 1).times(that.wij(i, 2))));
        case 5:
            return that.pik(i, 4)
                   .plus(Complex.times(1 / 5, that.wij(i, 2).times(that.wij(i, 2))));
        default:
            should.fail("invalid k:" + k);
        }
    };

#endif
