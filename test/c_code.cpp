/* Small C file creating an array to demo C -> Python data passing 
 * 
 * Author: Gael Varoquaux
 * License: BSD
 */

#include <stdlib.h>


void LiberaMat(double **Mat, int i){ 
    int k;
    for(k=0; k<=i; k++)
        free(Mat[i]);
    free(Mat);
}


double **AllocMat_(int nFilas, int nColumnas){
    double **Mat;

    if((Mat = (double**) calloc(nFilas, sizeof(double *))) == NULL)//porq "Doub *"?
        return NULL;

    int i;
    for(i=0; i<nFilas; i++)
        if((Mat[i] = (double *) calloc(nColumnas, sizeof(double)))==NULL) {
            LiberaMat(Mat,i);
            return NULL;
        }   
    return Mat;
}


double **AllocMat(int nRows, int nCols){
    double **Mat;
    // dsps probar con calloc!
    if((Mat = (double **) calloc(nRows,sizeof(double *))) == NULL)//porq "double *"? 
        return NULL;

    //--- alocamos el punto de origen
    if((Mat[0] = (double *) calloc( nCols*nRows,sizeof( double ) ))==NULL){
        LiberaMat(Mat,0);
        return NULL;
    }   

    //--- alojamos los demas elementos
    int i;
    for(i=1; i<nRows; i++)
        Mat[i] = Mat[0] + i*nCols; //(double *) malloc(nCols*sizeof(double));
    return Mat;
}


double *compute(int size){
    double* array;
    //array = malloc(sizeof(double)*size); // OK with C
    array = (double*) malloc(sizeof(double)*size);  // ok with c++
    int i;
    for (i=0; i<size; i++){
        array[i] = i;
    }
    return array;
}


double **compute_2d(int nx, int ny){
    double** array;
    array = AllocMat(nx, ny);
    int i, j;
    for (i=0; i<nx; i++){
        for (j=0; j<ny; j++){
            array[i][j] = 22. + 1.*(i+j);
        }
    }
    return array;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++
class array2d{
    private:
        double* m;

    public:
        int nrows, ncols;
        array2d(int nr, int nc): nrows(nr), ncols(nc) { 
            m = (double*) malloc(nr*nc*sizeof(double));
        }
        double* operator[](const int i){ 
            return &m[i*ncols]; 
        } //i'th element (must be a pointer)
};
/*
inline T & NRvector<T>::operator[](const int i){ //subscripting
    return &m[i*ncols];
}*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++


double *compute_2d_ii(int nx, int ny){
    array2d array(nx, ny);
    //double *array;
    //array = (double*) malloc(nx*ny*sizeof(double));
    int i, j;
    for (i=0; i<nx; i++){
        for (j=0; j<ny; j++){
            //array[i*ny+j] = 1.*(i+j);
            array[i][j] = 1.*(i+j);
        }
    }
    return &(array[0][0]); //&array[0];
}

/*
double **compute_2d(int nx, int ny){
    double** array;
    array = 
    int i;
    for (i=0; i<size; i++){
        array[i] = i;
    }
    return array;
}
*/

