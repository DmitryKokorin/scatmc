#ifndef _COMMON_H_
#define _COMMON_H_


#define DOUBLE_PRECISION 1

#if defined DOUBLE_PRECISION
#define Float double
#else
#define Float float
#endif



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



#endif /* _COMMON_H_ */
