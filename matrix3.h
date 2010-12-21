#ifndef MATRIX3_Z7HSXRKZ

#define MATRIX3_Z7HSXRKZ

#include "vector3.h"

class Matrix3 {
public:

	Matrix3 ()
	{
    	data[0][0] = 0.;
    	data[0][1] = 0.;
    	data[0][2] = 0.;
    	data[1][0] = 0.;
    	data[1][1] = 0.;
    	data[1][2] = 0.;
    	data[2][0] = 0.;
    	data[2][1] = 0.;
    	data[2][2] = 0.;
	}


	inline Float&   operator() (const int& row, const int& col)
	{
		return data[row][col];
	}

	Float    det() const;
	Vector3  operator*(const Vector3& rhv) const;
	Matrix3  operator*(Matrix3& rhv) const;

	friend inline Matrix3  invert(const Matrix3& mtx);


private:
	Float data[3][3];
};

#include "matrix3_inl.h"

#endif /* end of include guard: MATRIX3_Z7HSXRKZ */
