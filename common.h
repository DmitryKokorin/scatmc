#ifndef _COMMON_H_
#define _COMMON_H_


#define DOUBLE_PRECISION 1

#if defined DOUBLE_PRECISION
#define Float double
#else
#define Float float
#endif

extern const Float kMachineEpsilon;

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


template < typename T >
T **allocate2dArray( int rows, int cols)
{
    T **ppi = new T*[rows];
    T *cur_ptr = new T [rows * cols];

    for ( int i = 0; i < rows; ++i) {

        *(ppi + i) = cur_ptr;
         cur_ptr += cols;
    }

    return ppi;
}


template < typename T >
void free2dArray(T** array)
{
    delete [] *array;
    delete [] array;
}

int solveQuadric(const Float a, const Float b, const Float c, Float& x1, Float& x2);



#endif /* _COMMON_H_ */
