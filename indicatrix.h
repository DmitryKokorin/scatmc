#ifndef _INDICATRIX_H_
#define _INDICATRIX_H_

#include "direction.h"
#include "vector3.h"


class Indicatrix
{
public:

    Indicatrix(const Direction& d_i);
    virtual ~Indicatrix();

	Float operator()(const Direction& d_s);


private:
	
	Direction d_i;
	Vector3   k_i;
	Vector3   ee_i;
	Float     een_i;
	Float     factor1;

	static const Float kCalculationEpsilon;
};


#endif /* _INDICATRIX_H_ */
