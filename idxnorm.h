#ifndef _IDXNORM_H_

#define _IDXNORM_H_


namespace idxnorm {

template <class T>
class phiFunctor
{
public:

    phiFunctor(Indicatrix<T, Optics::OBeam>& indO_, Indicatrix<T, Optics::EBeam>& indE_,
                const Vector3& nn_, const Angle& a_i_,
                const Float sint_s_, const Float cost_s_) :
        indO(indO_),
        indE(indE_),
        nn(nn_),
        a_i(a_i_),
        sint_s(sint_s_),
        cost_s(cost_s_)
    {}

    inline Float operator()(const Float phi)
    {  
        Vector3 s_s = Vector3(sint_s*cos(phi), sint_s*sin(phi), cost_s);

		Angle a_s   = Angle(s_s, nn);

		return  sint_s * (  indO(s_s) + indE(s_s) );
    }

protected:

    Indicatrix<T, Optics::OBeam>& indO;
    Indicatrix<T, Optics::EBeam>& indE;

    const Vector3& nn;
    const Angle& a_i;

    Float sint_s;
    Float cost_s;
};


template <class T>
class thetaFunctor
{
public:

    thetaFunctor(Indicatrix<T, Optics::OBeam>& indO_, Indicatrix<T, Optics::EBeam>& indE_,
                    const Vector3& nn_, const Angle& a_i_) :
        indO(indO_),
        indE(indE_),
        nn(nn_),
        a_i(a_i_)
    {}

    inline Float operator()(const Float theta)
    {  
        Float cost_s = cos(theta);
        Float sint_s = sin(theta);
     
        phiFunctor<T> functor = phiFunctor<T>(indO, indE, nn, a_i, sint_s, cost_s);
        Adapt s(1.0e-10);
        return s.integrate(functor, 0., 2*M_PI);
    }

protected:

    Indicatrix<T, Optics::OBeam>& indO;
    Indicatrix<T, Optics::EBeam>& indE;

    const Vector3& nn;
    const Angle& a_i;
};
 
} //namespace




template <class T>
void createIndicatrixNorm(LinearInterpolation& li, const int kPoints = 1000)
{

    li.resize(kPoints);
    li.setRange(0., 0.5*M_PI);

	const Float kResolution = li.resolution();

    #pragma omp parallel for
	for (int i = 1; i < kPoints; ++i) {

		Float theta_i     = i*kResolution;
		Angle   a_i = Angle(theta_i);
		Vector3 s_i = createSomeDeviantVector(Optics::director, a_i).normalize();

		//create coordinate system

		Vector3 v2 = crossProduct(s_i, Optics::director).normalize();
		Vector3 v3 = crossProduct(v2, s_i);

		Matrix3 mtx = invert(createTransformMatrix(s_i, v2, v3));
		Vector3 nn = mtx*Optics::director;

        Indicatrix<T, Optics::OBeam> indO = Indicatrix<T, Optics::OBeam>(Vector3(1., 0., 0.), nn);
		Indicatrix<T, Optics::EBeam> indE = Indicatrix<T, Optics::EBeam>(Vector3(1., 0., 0.), nn);

        Adapt s(1.0e-10);
        idxnorm::thetaFunctor<T> functor(indO, indE, nn, a_i);
        Float integral = s.integrate(functor, 0., M_PI);

		li[i] = integral;

		fprintf(stderr, "%d\n", i);
	}

	li[0] = li[1];
}




#endif /* end of include guard: _IDXNORM_H_ */

