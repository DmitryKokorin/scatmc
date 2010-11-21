#ifndef MATRIX3_Z7HSXRKZ

#define MATRIX3_Z7HSXRKZ

#include "vector3.h"

class Matrix3 {
public:
	Matrix3 ();

	inline Float&   operator() (const int& row, const int& col);

	inline Float    det() const;
	inline Vector3  operator*(const Vector3& rhv) const;

	friend inline Matrix3  invert(const Matrix3& mtx);


private:
	Float data[3][3];
};

#endif /* end of include guard: MATRIX3_Z7HSXRKZ */
