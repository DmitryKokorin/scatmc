#ifndef _EXTLENGTH_H_ 
#define _EXTLENGTH_H_

#include <vector>

#include "common.h"

class Angle;


class ExtLength 
{
public:
    ExtLength(const int kThetaIterations = 1000, const int kPhiIterations = 1000);
    virtual ~ExtLength();

	Float operator()(const Angle& d) const;

	static const int    kPoints = 400;
	static const Float  kResolution;

	//private:
	Float lengths[kPoints]; //plain array to have this data in a stack and so in a cache 
};


#endif /* _EXTLENGTH_H_ */
