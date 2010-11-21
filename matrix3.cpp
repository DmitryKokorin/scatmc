#include <memory.h>
#include <cmath>
#include "matrix3.h"


Matrix3::Matrix3()
{
	memset(&data, 0, sizeof(data));
}

Float& Matrix3::operator()(const int& row, const int& col)
{
	return data[row][col];
}

Float Matrix3::det() const
{
	return +data[0][0]*(data[1][1]*data[2][2]-data[2][1]*data[1][2])
		   -data[0][1]*(data[1][0]*data[2][2]-data[1][2]*data[2][0])
		   +data[0][2]*(data[1][0]*data[2][1]-data[1][1]*data[2][0]);
}

Vector3 Matrix3::operator*(const Vector3& rhv) const
{
	Float x = data[0][0]*rhv.x() + data[0][1]*rhv.y() + data[0][2]*rhv.z();
	Float y = data[1][0]*rhv.x() + data[1][1]*rhv.y() + data[1][2]*rhv.z();
	Float z = data[2][0]*rhv.x() + data[2][1]*rhv.y() + data[2][2]*rhv.z();

	return Vector3(x, y, z);
}


Matrix3 invert(const Matrix3& mtx) 
{
	Float id = 1. / mtx.det();

	if (isnan(id))
		return Matrix3(); //singular matrix, return something


	Matrix3 res;

	res(0,0) =  (mtx.data[1][1]*mtx.data[2][2]-mtx.data[2][1]*mtx.data[1][2])*id;
	res(1,0) = -(mtx.data[0][1]*mtx.data[2][2]-mtx.data[0][2]*mtx.data[2][1])*id;
	res(2,0) =  (mtx.data[0][1]*mtx.data[1][2]-mtx.data[0][2]*mtx.data[1][1])*id;
	res(0,1) = -(mtx.data[1][0]*mtx.data[2][2]-mtx.data[1][2]*mtx.data[2][0])*id;
	res(1,1) =  (mtx.data[0][0]*mtx.data[2][2]-mtx.data[0][2]*mtx.data[2][0])*id;
	res(2,1) = -(mtx.data[0][0]*mtx.data[1][2]-mtx.data[1][0]*mtx.data[0][2])*id;
	res(0,2) =  (mtx.data[1][0]*mtx.data[2][1]-mtx.data[2][0]*mtx.data[1][1])*id;
	res(1,2) = -(mtx.data[0][0]*mtx.data[2][1]-mtx.data[2][0]*mtx.data[0][1])*id;
	res(2,2) =  (mtx.data[0][0]*mtx.data[1][1]-mtx.data[1][0]*mtx.data[0][1])*id;

	return res;
}
