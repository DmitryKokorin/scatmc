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
	x1(0), x2(0), y1(0), y2(0), width(0), height(0), square(0)
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
