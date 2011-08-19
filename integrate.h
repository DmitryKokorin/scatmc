#ifndef _INTEGRATE_H_
#define _INTEGRATE_H_

#include <cstdio>
#include "mathcompat.h"
#include "common.h"


template <class T>
class Trapzd
{
public:

    Trapzd(T& func_, const Float l_, const Float r_) :
        n(0),
        func(func_),
        l(l_),
        r(r_),
        s(0.)
    {}

    Float next()
    {
        Float x, tnm, sum, del;
        int it,j;
        ++n;

        if (1 == n) {

            return s = 0.5*(r - l)*(func(l) + func(r));
        }
        else {

            for (it = 1, j = 1; j < n-1; ++j)
                it <<= 1;

            tnm = it;
            del = (r - l) / tnm;

            x = l + 0.5*del;

            for (sum = 0., j = 0; j < it; ++j, x += del)
                sum += func(x);

            s = 0.5*(s + del*sum);

            return s;
        }
    }

protected:

    int n;

    T&    func;

    Float l, r;
    Float s;
};


template<class T>
Float qsimp(T &func, const Float l, const Float r, const Float eps)
{
    const int JMAX = 20;
    Float s, st, ost = 0.0, os = 0.0;
    Trapzd<T> t(func, l, r);

    for (int j = 0; j < JMAX; ++j) {

        st = t.next();
        s  = (4.0*st - ost)/3.0;

        if (j > 5) {
            if (fabs(s-os) < eps*fabs(os) || (s == 0.0 && os == 0.0)) {

                return s;
            }
        }

        os  = s;
        ost = st;
    }

    fprintf(stderr, "Too many steps in routine qsimp\n");
    return s;
}




class Adapt
{
public:

	Adapt(const Float tol);

	template <class T>
	Float integrate(T &func, const Float a, const Float b);

private:

   	template <class T>
	Float adaptlob(T &func, const Float a, const Float b, const Float fa, const Float fb, const Float is);

    static const Float alpha, beta, x1, x2, x3, x[12];

	Float TOL,toler;
	
	bool terminate;
	bool out_of_tolerance;
};


template <class T>
Float Adapt::integrate(T &func, const Float a, const Float b)
{
	Float y[13];

	Float m = 0.5*(a + b);
	Float h = 0.5*(b - a);

	Float fa = y[0]  = func(a);
	Float fb = y[12] = func(b);

	for (int i = 1; i < 12; ++i)
		y[i] = func(m + x[i]*h);

	Float i2 = (h/6.0)*(y[0] + y[12] + 5.0*(y[4] + y[8]));

	Float i1 = (h/1470.0)*(77.0  *(y[0]+y[12])+
	                       432.0 *(y[2]+y[10])+
		                   625.0 *(y[4]+y[8]) +
    		               672.0 *y[6]);

	Float is = h*(0.0158271919734802*(y[0]+y[12]) + 
	              0.0942738402188500*(y[1]+y[11]) + 
	              0.155071987336585 *(y[2]+y[10]) +
    		      0.188821573960182 *(y[3]+y[9])  +
           		  0.199773405226859 *(y[4]+y[8])  + 
        		  0.224926465333340 *(y[5]+y[7])  +
    		      0.242611071901408 *y[6]);

	Float erri1 = fabs(i1-is);
	Float erri2 = fabs(i2-is);

	Float r = (erri2 != 0.0) ? erri1/erri2 : 1.0;
	toler = (r > 0.0 && r < 1.0) ? TOL/r : TOL;

	if (0.0 == is)
		is = b - a;

	is = fabs(is);

	return adaptlob(func, a, b, fa, fb, is);
}

template <class T>
Float Adapt::adaptlob(T &func, const Float a, const Float b,
                        const Float fa, const Float fb, const Float is)
{
	Float m = 0.5*(a + b);
	Float h = 0.5*(b - a);

	Float mll = m - alpha*h;
	Float ml  = m - beta*h;
	Float mr  = m + beta*h;
	Float mrr = m + alpha*h;

	Float fmll = func(mll);
	Float fml  = func(ml);
	Float fm   = func(m);
	Float fmr  = func(mr);
	Float fmrr = func(mrr);

	Float i2 = h/6.0*(fa + fb + 5.0*(fml + fmr));

	Float i1 = h/1470.0*(77.0*(fa   + fb)  +
	                    432.0*(fmll + fmrr)+
	                    625.0*(fml  + fmr) +
	                    672.0*fm);

	if (fabs(i1-i2) <= toler*is || mll <= a || b <= mrr) {

		if ((mll <= a || b <= mrr) && terminate) {

			out_of_tolerance = true;
			terminate = false;
		}

		return i1;
	}
	else
		return  adaptlob(func, a,   mll, fa,   fmll, is) +
    			adaptlob(func, mll, ml,  fmll, fml,  is) +
			    adaptlob(func, ml,  m,   fml,  fm,   is) +
			    adaptlob(func, m,   mr,  fm,   fmr,  is) +
			    adaptlob(func, mr,  mrr, fmr,  fmrr, is) +
			    adaptlob(func, mrr, b,   fmrr, fb,   is);
}



#endif /* end of include guard: _INTEGRATE_H_ */
