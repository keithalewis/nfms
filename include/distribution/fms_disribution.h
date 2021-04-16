// fms_distribution.h - Distributions and derivatives
#pragma once

namespace fms {

	template<class X = double, class S = double>
	class distribution {
	public:
		using xtype = typename X;
		using stype = typename S;

		// d^n/dx^n P_s(X <= x)
		X cdf(const X& x, const S& s, size_t n) const
		{
			return cdf_(x, s, n);
		}
	private:
		virtual X cdf_(const X&, const S&, size_t) = 0;
	}

} // namespace fms
