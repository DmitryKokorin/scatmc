#ifndef _NODE_H_
#define _NODE_H_



class GreedRect
{
public:
	GreedRect();
	GreedRect(const int& x1_, const int& y1_, const int& x2_, const int& y2_) :
		x1(x1_),
		y1(y1_),
		x2(x2_),
		y2(y2_),
		midx((x1+x2)/2),
		midy((y1+y2)/2),
		width(x2-x1),
		height(y2-y1),
		square(width*height)
	{
	}

	GreedRect topHalf()    const { return GreedRect(x1, y1, x2, midy); }
	GreedRect bottomHalf() const { return GreedRect(x1, midy, x2, y2); }
	GreedRect leftHalf()   const { return GreedRect(x1, y1, midx, y2); }
	GreedRect rightHalf()  const { return GreedRect(midx, y1, x2, y2); }

	bool canSplitX()   const    { return x1 != midx; }
	bool canSplitY()   const    { return y1 != midy; }

	int x1, y1, x2, y2;
	int midx, midy;
	int width, height;
	int square;
};



class Node
{
public:

	Node(const GreedRect& rect_);
	~Node();

	bool isLeaf();

	bool splitX();
    bool splitY(); 

	GreedRect rect;


	Node* pParent;
	Node* pChild1;
	Node* pChild2;

	
private:	
	//disable copying
	Node(const Node&);
	Node& operator=(const Node&);
};


#endif /* _NODE_H_ */
