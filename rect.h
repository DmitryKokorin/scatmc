#ifndef _RECT_H_
#define _RECT_H_

#include <vector>
#include "common.h"

struct Knot
{
	Knot(const Float x_, const Float y_, const Float val_) :
		x(x_),
		y(y_),
		val(val_)
	{};

	Knot() :
		x(0.),
		y(0.),
		val(0.)
	{};

	Float x;
	Float y;
	Float val;
};

struct Rect
{
	Rect(	const int tl, const int tr,
			const int bl, const int br);
	
	int	tl,	tr,
		bl,	br;

	Float	x1, x2, y1, y2;

	Float	width,
			height;
			
	Float	square;

	
	inline 	Float	integral();

	static std::vector<Knot>	*s_knots;
};

#endif /* _RECT_H_ */
