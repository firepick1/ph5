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
PHFeed<T>::PHFeed(PH5Curve<T> &ph5, T vMax, T tvMax, T vIn, T vOut, T vCruise) 
	:ph(ph5), vMax(vMax), tvMax(tvMax), vIn(vIn), vCruise(vCruise), vOut(vOut) {
	S = ph.s(1);
	ASSERT(vMax > 0);
	ASSERT(tvMax > 0);
	ASSERT(0 <= vIn && vIn <= vMax);
	if (vCruise == 0) {
		vCruise = vMax;
	}
	ASSERT(0 <= vCruise && vCruise <= vMax);
	ASSERT(0 <= vOut && vOut <= vMax);
	iterations = 10;
	epsilon = 0.0000001;
	T sMax = vMax * tvMax/2.0; // distance required to reach max velocity
	if (vIn == vCruise && vCruise == vOut) { // constant velocity
		uc = FEEDUSE_A;
		sAccel = 0;
		tAccel = 0;
		sDecel = 0;
		tDecel = 0;
	} else if (vIn != vCruise && vCruise == vOut) {	// acceleration only
		if (sMax > S) {
			uc = FEEDUSE_B1;
			T sRatio = sqrt(S/sMax);
			sAccel = S;
			tAccel = tvMax * sRatio;
			vCruise = vCruise * sRatio;
		} else {
			uc = FEEDUSE_B2;
			sAccel = sMax;
			tAccel = tvMax;
		}
		sDecel = 0;
		tDecel = 0;
	} else if (vIn == vCruise && vCruise != vOut) {	// deceleration only
		if (sMax > S) {
			uc = FEEDUSE_C1;
			T sRatio = sqrt(S/sMax);
			vCruise = vCruise * sRatio;
			sDecel = S;
			tDecel = tvMax * sRatio;
		} else {
			uc = FEEDUSE_C2;
			sDecel = sMax;
			tDecel = tvMax;
		}
		sAccel = 0;
		tAccel = 0;
	} else {	// accelerate, cruise, decelerate
		T S2 = S/2.0;
		ASSERTEQUAL(vIn, vOut);
		if (sMax > S2) {
			uc = FEEDUSE_D1;
			sAccel = S2;
			T sRatio = sqrt(S2/sMax);
			vCruise *= sRatio;
			tAccel = tvMax * sRatio;
		} else {
			uc = FEEDUSE_D2;
			sAccel = sMax;
			tAccel = tvMax;
		}
		sDecel = sAccel;
		tDecel = tAccel;
	}
	sCruise = S - sAccel - sDecel;
	tCruise = sCruise/vCruise;
	tS = tAccel + tCruise + tDecel;
	tauCruise = tAccel/tS;
	tauDecel = 1 - tDecel/tS;

	for (int16_t i=0; i<=6; i++) {
		Faccel[i] = Fk(vIn, vCruise, i);
		Fcruise[i] = Fk(vCruise, vCruise, i);
		Fdecel[i] = Fk(vCruise, vOut, i);
	}
}

template <class T>
T PHFeed<T>::Ekt(T Ekprev, T tau) {
	T Ekr = Ekprev;
	T Ftau = Ft(tau);
	T dE = 0;

	for (int16_t iteration=0; iteration<iterations; iteration++) {
		dE = (Ftau - ph.s(Ekr))/ph.sigma(Ekr);
		Ekr = Ekr + dE;
		//cout << "tau:" << tau << " Ftau:" << Ftau << " Ekr:" << Ekr << " dE:" << dE << endl;
		if (Ekr < 0) { Ekr = 0; }
		if (Ekr > 1) { Ekr = 1; }
		if (-epsilon < dE && dE < epsilon) {
			return Ekr;
		}
	}
	LOGDEBUG1("Ekt() exceeded iterations dE:%lf", (double) dE);
	return Ekr;
}

template <class T>
T PHFeed<T>::Ft(T tau) {
	T sum = 0;
	if (tau < tauCruise) {			// accelerating
		T t = tau ? (tau*tS)/tAccel : 0;
		for (int16_t k=0; k<=6; k++) {
			sum += Faccel[k] * Bernstein6(k, t);
		}
		return sum*tAccel/6.0;
	} else if (tau < tauDecel) {	// cruising
		T t = (tau*tS-tAccel)/tCruise;
		for (int16_t k=0; k<=6; k++) {
			sum += Fcruise[k] * Bernstein6(k, t);
		}
		return sum*tCruise/6.0 + sAccel;
	} else {						// decelerating
		T t = tau==1 ? 1 : (tau*tS-tAccel-tCruise)/tDecel;
		for (int16_t k=0; k<=6; k++) {
			sum += Fdecel[k] * Bernstein6(k, t);
		}
		return sum*tDecel/6.0 + sAccel + sCruise;
	}
}

template <class T>
T PHFeed<T>::Fk(T vIn, T vOut, int16_t k) {
	T sum = 0;
	for (int16_t j=0; j<k; j++) {
		sum += j<3 ? vIn : vOut;
	}
	return sum;
}

template class PHFeed<PH5TYPE>;
