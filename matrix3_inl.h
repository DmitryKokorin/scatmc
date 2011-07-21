#include "mathcompat.h"
#include "matrix3.h"



inline Float Matrix3::det() const
{
	return +data[0][0]*(data[1][1]*data[2][2]-data[2][1]*data[1][2])
		   -data[0][1]*(data[1][0]*data[2][2]-data[1][2]*data[2][0])
		   +data[0][2]*(data[1][0]*data[2][1]-data[1][1]*data[2][0]);
}

inline Vector3 Matrix3::operator*(const Vector3& rhv) const
{
	Float x = data[0][0]*rhv.x() + data[0][1]*rhv.y() + data[0][2]*rhv.z();
	Float y = data[1][0]*rhv.x() + data[1][1]*rhv.y() + data[1][2]*rhv.z();
	Float z = data[2][0]*rhv.x() + data[2][1]*rhv.y() + data[2][2]*rhv.z();

	return Vector3(x, y, z);
}


inline Matrix3 invert(const Matrix3& mtx)
{
	Float id = 1. / mtx.det();

	if (isnan(id))
		return Matrix3(); //singular matrix, return something


	Matrix3 res;

	res(0,0) =  (mtx.data[2][2]*mtx.data[1][1] - mtx.data[2][1]*mtx.data[1][2])*id;
	res(0,1) = -(mtx.data[2][2]*mtx.data[0][1] - mtx.data[2][1]*mtx.data[0][2])*id;
	res(0,2) =  (mtx.data[1][2]*mtx.data[0][1] - mtx.data[1][1]*mtx.data[0][2])*id;

	res(1,0) = -(mtx.data[2][2]*mtx.data[1][0] - mtx.data[2][0]*mtx.data[1][2])*id;
	res(1,1) =  (mtx.data[2][2]*mtx.data[0][0] - mtx.data[2][0]*mtx.data[0][2])*id;
	res(1,2) = -(mtx.data[1][2]*mtx.data[0][0] - mtx.data[1][0]*mtx.data[0][2])*id;

	res(2,0) =  (mtx.data[2][1]*mtx.data[1][0] - mtx.data[2][0]*mtx.data[1][1])*id;
	res(2,1) = -(mtx.data[2][1]*mtx.data[0][0] - mtx.data[2][0]*mtx.data[0][1])*id;
	res(2,2) =  (mtx.data[1][1]*mtx.data[0][0] - mtx.data[1][0]*mtx.data[0][1])*id;


	return res;
}

inline Matrix3 Matrix3::operator*(Matrix3& rhv) const
{
	Matrix3 res;

	for (int j = 0; j < 3; ++j)
		for (int i = 0; i < 3; ++i)
			for (int k = 0; k < 3; ++k)
				res(j,i) += data[j][k]* rhv(k,i);

	return res;
}
