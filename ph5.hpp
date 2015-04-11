#ifndef PH5_HPP
#define PH5_HPP

#include <vector>
#include <map>
#ifdef _MSC_VER
#include "winjunk.hpp"
#else
#define CLASS_DECLSPEC
#endif
#include "jansson.h"
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

CLASS_DECLSPEC typedef map<string, const char *> ArgMap;
CLASS_DECLSPEC extern ArgMap emptyMap;

template <class T>
class Complex {
	private:
		T re;
		T im;

	public:
		Complex(T re=0, T im=0);
	public:
		inline T Re() { return this->re; }
	public:
		inline T Im() { return this->im; }
	public: 
		T modulus() { return ::sqrt(re*re + im*im); }
	public:
		inline Complex<T> operator+(Complex<T> &that) { 
			return Complex(re+that.re, im+that.im);
		}
	public:
		inline Complex<T> operator+(T k) { 
			return Complex(re+k, im);
		}
	public:
		inline Complex<T> operator-(Complex<T> &that) { 
			return Complex(re-that.re, im-that.im);
		}
	public:
		inline Complex<T> operator-(T k) { 
			return Complex(re-k, im);
		}
	public:
		inline bool operator==(Complex<T> that) {
			return re == that.re && im == that.im;
		}
	public:
		inline Complex<T> operator*(Complex<T> &that) {
			return Complex(re*that.re-im*that.im,re*that.im+im*that.re);
		}
	public:
		inline Complex<T> operator*(T that) {
			return Complex(re*that,im*that);
		}
	public:
		inline Complex<T> conj() {
			return Complex(re, -im);
		}
	public:
		inline Complex<T> recip() {
			T denom = re*re + im*im;
			return Complex(re/denom, -im/denom);
		}
	public:
		inline Complex<T> operator/ (Complex<T> &c2) {
			T denom = c2.re*c2.re + c2.im*c2.im;
			return Complex(
				(re*c2.re + im*c2.im)/denom,
				(im*c2.re - re*c2.im)/denom
			);
		}
	public:
		inline Complex<T> operator/ (T k) {
			return Complex(re/k, im/k);
		}
	public:
		inline void add(Complex<T> that) {
			re += that.re;
			im += that.im;
		}
	public:
		inline Complex<T> sqrt() {
			T m = modulus();
			T p = ::sqrt((m+re)/2);
			T q = ::sqrt((m-re)/2);
			if (im >= 0) {
				return Complex(p,q);
			} else {
				return Complex(p,-q);
			}
		}
	public:
		string stringify(int nPlaces=0);
};

namespace ph5 {
template<class T>
inline Complex<T> operator*(T k, Complex<T> c) {
	return Complex<T>(k*c.re, k*c.im);
}
}

typedef class CLASS_DECLSPEC ASDF {
    public:

} ASDF;

} // namespace ph5

#endif
