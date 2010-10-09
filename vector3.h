#ifndef _VECTOR3_H_
#define _VECTOR3_H_

#include "common.h"
//#include "direction.h"


class Vector3
{
public:
    Vector3();
	Vector3(const Vector3& rhv);
	Vector3(const Float x, const Float y, const Float z);

	inline Vector3& operator= (const Vector3& rhv);
	inline Vector3& operator+=(const Vector3& rhv);
	inline Vector3  operator+ (const Vector3& rhv) const;
	inline Vector3& operator-=(const Vector3& rhv);
	inline Vector3  operator- (const Vector3& rhv) const;

	inline Vector3& operator*=(const Float& rhv);
	inline Vector3  operator* (const Float& rhv) const;
   	inline Vector3& operator/=(const Float& rhv);
	inline Vector3  operator/ (const Float& rhv) const;

	inline Float    operator*(const Vector3& rhv) const;
   
	inline Float x() const; 
	inline Float y() const;
	inline Float z() const;

	inline Float norm() const;
	inline Vector3& normalize();

private:
	Float x_, y_, z_;
};

inline Vector3 operator*(const Float& lhv, const Vector3& rhv);
inline Vector3 crossProduct(const Vector3& lhv, const Vector3& rhv);


#include "vector3_inl.h"


#endif /* _VECTOR3_H_ */
