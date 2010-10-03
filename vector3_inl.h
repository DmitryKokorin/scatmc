#ifndef _VECTOR3_INL_H_
#define _VECTOR3_INL_H_

#include <cmath>


Float Vector3::x() const
{
	return x_;
}

Float Vector3::y() const
{
	return y_;
}

Float Vector3::z() const
{
	return z_;
}

Vector3& Vector3::operator=(const Vector3& rhv)
{
	if (this != &rhv) {
		x_ = rhv.x_;
		y_ = rhv.y_;
		z_ = rhv.z_;
	}
	
	return *this;
}

Vector3& Vector3::operator+=(const Vector3& rhv)
{
	x_ += rhv.x_;
	y_ += rhv.y_;
	z_ += rhv.z_;
	
	return *this;
}

Vector3 Vector3::operator+(const Vector3& rhv) const
{
	return Vector3(*this)+=rhv;
}

Vector3& Vector3::operator-=(const Vector3& rhv)
{
	x_ -= rhv.x_;
	y_ -= rhv.y_;
	z_ -= rhv.z_;
	
	return *this;
}

Vector3 Vector3::operator-(const Vector3& rhv) const
{
	return Vector3(*this)-=rhv;
}


Vector3& Vector3::operator*=(const Float& rhv)
{
	x_ *= rhv;
	y_ *= rhv;
	z_ *= rhv;
	
	return *this;
}

Vector3 Vector3::operator*(const Float& rhv) const
{
	return Vector3(*this)*=rhv;
}


Vector3& Vector3::operator/=(const Float& rhv)
{
	Float irhv = 1. / rhv;

	return *this *= irhv;
}

Vector3 Vector3::operator/(const Float& rhv) const
{
	return Vector3(*this)/=rhv;
}

Float Vector3::norm() const
{
	return sqrt(x_*x_ + y_*y_ + z_*z_);
}

Vector3& Vector3::normalize()
{
	Float n = norm();

	if (0.0 == n) {    //TODO: use some epsilon

		return *this;
	}

	return *this /= n;
}

Vector3 operator*(const Float& lhv, const Vector3& rhv)
{
	return rhv*lhv;
}

Vector3 crossProduct(const Vector3& lhv, const Vector3& rhv)
{
	return Vector3((lhv.y()*rhv.z() - lhv.z()*rhv.y()),
				  (-lhv.x()*rhv.z() + lhv.z()*rhv.x()),
				   (lhv.x()*rhv.y() - lhv.x()*rhv.y()));
}

#endif /* _VECTOR3_INL_H_ */
