#ifndef _INDICATRIX_H_
#define _INDICATRIX_H_

#include "vector3.h"
#include "angle.h"


class Indicatrix
{
public:

    Indicatrix(const Vector3& s_i, const Vector3& n_i);
    virtual ~Indicatrix();

	Float operator()(const Vector3& s_s);


private:
	
	Vector3   s_i;
	Vector3   k_i;
	Vector3   n_i; //director in k_i based coordinate system
	Angle     a_i; //angle between k_i and n vectors 

	Vector3   ee_i;
	Float     een_i;
	Float     factor1;

	static const Float kCalculationEpsilon;
};


#endif /* _INDICATRIX_H_ */
