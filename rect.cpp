#include <cstdio>
#include "rect.h"

std::vector<Knot> *Rect::s_knots = NULL;

Rect::Rect(const int tl_,
		   const int tr_,
		   const int bl_,
		   const int br_) :
	tl(tl_),
	tr(tr_),
	bl(bl_),
	br(br_),
	x1(0), x2(0), y1(0), y2(0), width(0), height(0), square(0),
	val(0)
{
	std::vector<Knot>& knots = *s_knots;
	x1 = knots[tl].x;
	x2 = knots[tr].x;
	y1 = knots[tl].y;
	y2 = knots[bl].y;

	width  = x2 - x1;
	height = y2 - y1;

	square = width*height;
}

void Rect::choosePointInRect(Float& x, Float& y, const Float randX, const Float /*randY*/)
{
	std::vector<Knot>& knots = *s_knots;

	Float b1 = knots[tl].val;
	Float b2 = knots[tr].val - b1;
	Float b3 = knots[bl].val - b1;
	Float b4 = b1 - knots[tr].val - knots[bl].val + knots[br].val;


	//fprintf(stderr, "tl=%f\ttr=%f\tbl=%f\tbr=%f\n", knots[tl].val, knots[tr].val, knots[bl].val, knots[br].val);
	//fprintf(stderr, "b1=%f\tb2=%f\tb3=%f\tb4=%f\n", b1, b2, b3, b4);

	int roots;

	{
		Float x1, x2;
		roots = solveQuadric(b2 + 0.5*b4, 0.5*(b1 + 0.5*b3), -randX*(b2+0.5*b4+0.5*b1+0.25*b3), x1, x2);

		if (roots == 1)
			x = x1;
		else if (roots == 2) {

			if (x1 >= 0. && x1 < 1.)
				x = x1;
			else if (x2 >= 0. && x2 < 1.)
				x = x2;
			else
				fprintf(stderr, "x out of range, %f\t%f\n", x1, x2);
		}
	}

	{
		Float y1, y2;
		roots = solveQuadric(b3 + b4*x, 0.5*(b1 + b2*x), /*-randY*(b1+b2*x + 0.5*(b3+b4*x)*/0., y1, y2);

		if (roots == 1)
			y = y1;
		else if (roots == 2) {

			if (y1 >= 0. && y1 < 1.)
				y = y1;
			else if (y2 >= 0. && y2 < 1.)
				y = y2;
			else
				fprintf(stderr, "y out of range, %f\t%f\n", y1, y2);
		}
	}

	x *= width;
	y *= height;

	x += knots[tl].x;
	y += knots[tl].y;
}


