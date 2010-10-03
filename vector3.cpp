#include "vector3.h"


Vector3::Vector3() :
	x_(0.0),
	y_(0.0),
	z_(0.0)
{
}

Vector3::Vector3(const Vector3& rhv) :
	x_(rhv.x_),
	y_(rhv.y_),
	z_(rhv.z_)
{
}

Vector3::Vector3(const Float x, const Float y, const Float z) :
	x_(x),
	y_(y),
	z_(z)
{
}
