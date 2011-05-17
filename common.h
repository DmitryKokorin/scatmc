#ifndef _COMMON_H_
#define _COMMON_H_


#define DOUBLE_PRECISION 1

#if defined DOUBLE_PRECISION
#define Float double
#else
#define Float float
#endif

extern const Float kMachineEpsilon;

//[y][x], x -- width, y -- height

template < typename T >
T **allocate2dArray(int height, int width)
{
    T **ppi = new T*[height];
    T *cur_ptr = new T [height * width];

    for ( int i = 0; i < height; ++i) {

        *(ppi + i) = cur_ptr;
         cur_ptr += width;
    }

    return ppi;
}

template < typename T >
void free2dArray(T** array)
{
    delete [] *array;
    delete [] array;
}


//[z][y][x], x -- width, y -- height, z -- depth

template < typename T >
T ***allocate3dArray(int depth, int height, int width)
{
    T *allElements = new T [depth * height * width];
    T ***array3D = new T** [depth];

    for (int i = 0; i < depth; ++i) {

        array3D[i] = new T*[height];

        for (int j = 0; j < height; ++j) {

            array3D[i][j] = allElements + (i * height * width) + (j * width);
        }
    }

    return array3D;
}

//TODO: fix this or allocation code
template < typename T >
void free3dArray(T*** array)
{
    delete [] **array;
    delete [] *array;
    delete [] array;
}



int solveQuadric(const Float a, const Float b, const Float c, Float& x1, Float& x2);



#endif /* _COMMON_H_ */
