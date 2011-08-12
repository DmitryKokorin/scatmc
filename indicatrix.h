#ifndef _INDICATRIX_H_
#define _INDICATRIX_H_

#include "vector3.h"
#include "angle.h"


class IndicatrixBase
{
public:

    IndicatrixBase(const Vector3& s_i, const Vector3& director_i);
    virtual ~IndicatrixBase() {}


protected:
	
	Vector3   s_i;
	Vector3   k_i;
	Vector3   director_i;   //director in k_i based coordinate system
	Angle     a_i;          //angle between k_i and n vectors 
	Vector3   e_i;

	Float     ei_n;
	Float     factor1;

	static const Float kCalculationEpsilon;
};

class IndicatrixEE : public IndicatrixBase
{
public:

    IndicatrixEE(const Vector3& s_i, const Vector3& director_i);
    virtual ~IndicatrixEE() {}

	Float operator()(const Vector3& s_s);
};

class IndicatrixEO : public IndicatrixBase
{
public:

    IndicatrixEO(const Vector3& s_i, const Vector3& director_i);

	Float operator()(const Vector3& s_s);
};

class IndicatrixOE : public IndicatrixBase
{
public:

    IndicatrixOE(const Vector3& s_i, const Vector3& director_i);

	Float operator()(const Vector3& s_s);
};

#ifdef EXPERIMENTAL

class IndicatrixHenyey : public IndicatrixBase
{
public:

    IndicatrixHenyey(const Vector3& s_i, const Vector3& director_i);

	Float operator()(const Vector3& s_s);
};

#endif


#endif /* _INDICATRIX_H_ */
