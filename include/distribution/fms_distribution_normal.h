// fms_distribution_normal.h - Standard normal distribution.
#pragma once
#include <cmath>
#include "fms_distribution.h"

namespace fms {
	
	template<class X, class S>
	class distribution_normal : public distribution<X,S> {
#ifndef  M_SQRT2
		static constexpr X M_SQRT2 = X(1.41421356237309504880);
#endif 
#ifndef  M_SQRT2PI
		static constexpr X M_SQRT2PI = X(2.50662827463100050240);
#endif 
	public:
		X cdf_(const X& x, const S& s, size_t n) const override
		{
			X x_ = x - s;

			if (n == 0) {
				return (1 + erf(x_ / X(M_SQRT2))) / 2;
			}

			X phi = exp(-x_ * x_ / X(2)) / X(M_SQRT2PI);

			if (n == 1) {
				return phi;
			}
														}
			// (d/dx)^n phi(x) = (-1)^n phi(x) H_n(x)
			return phi * H(n - 1, x_) * ((n&1) ? 1 : -1);
		}
	private:
		// Hermite polynomials
		constexpr X H(unsigned n, X x)
		{
			if (n == 0) {
				return X(1);
			}
			if (n == 1) {
				return x;
			}

			return x * H(n - 1, x) - X(n - 1) * H(n - 2, x);
		}
	};

} // namespace fms
