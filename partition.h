#ifndef _PARTITION_H_
#define _PARTITION_H_

struct Rect
{
	Rect(const Float& x1_, const Float& y1_, const Float& x2_, const Float& y2_) :
		x1(x1_), y1(y1_), x2(x2_), y2(y2_)

	Float x1, y1, x2, y2;
}

struct Node
{
	Node(const Rect& rect_) :
		pParent(NULL), pChild1(NULL), pChild2(NULL), rect(rect)
	{}

	~Node()
	{ 
		delete pChild1;
		delete pChild2;
	}

	Node* pParent;
	Node* pChild1;
	Node* pChild2;

	bool isLeaf() { return (NULL == pChild1) && (NULL == pChild2)}

	void splitX() 
	{ 
		pChild1 = new Node(Rect(x1, y1, 0.5 * (x1 + x2), y2));
		pChild2 = new Node(Rect(0.5 * (x1 + x2), y1, x2, y2));
		pChild1->pParent = this;
		pChild2->pParent = this;
	}

	void splitY() 
	{ 
		pChild1 = new Node(Rect(x1, y1, x2, 0.5*(y1+y2)));
		pChild2 = new Node(Rect(x1, 0.5*(y1+y2), x2, y2));
		pChild1->pParent = this;
		pChild2->pParent = this;
	}


	Rect rect;
};

enum TNodeStatus
{
	STAY_LEAF,
	SHOULD_SPLIT_X,
	SHOULD_SPLIT_Y,
};

int evaluateNodeStatus(const Node*& rect)
{
	return STAY_LEAF;
}

void createPartition(Node*& pRoot_)
{
	Node* pRoot = pRoot_;
	Node* pNode = pRoot;

	const Float epsilon = 1e-6;

	while () {
	}
}


#endif /* _PARTITION_H_ */
