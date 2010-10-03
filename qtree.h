#ifndef _QTREE_H_
#define _QTREE_H_


struct Rect
{
	Rect(Float x_, Float y_, Float width_, Float height_) :
		x(x_), y(y_), width(width_), height(height_);

	Float x, y, width, height;
};


class Node
{
public:
	Node();
	virtual ~Node();

	void subdivide();
	
protected:

	Rect  rect;
 
	Node* childTL;
	Node* childTR;
	Node* childBL;
	Node* childBR;
};


class QuadTree
{
public:
    QuadTree();
    virtual ~QuadTree();

};


#endif /* _QTREE_H_ */
