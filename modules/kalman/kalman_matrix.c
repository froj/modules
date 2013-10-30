#include "kalman_matrix.h"
#include <stdlib.h>
#include <stdio.h>

struct kalman_matrix{
	int cols;
	int rows;
	float* vals;
};

float* val(struct kalman_matrix* m, int p_col, int p_row);

void kalman_matrix_init(struct kalman_matrix** m,
                        int p_cols,
                        int p_rows,
                        float* p_vals){

    *m = (struct kalman_matrix*)malloc(sizeof(struct kalman_matrix));
    (*m)->cols = p_cols;
    (*m)->rows = p_rows;
    (*m)->vals = (float*)malloc(sizeof(float) * p_cols * p_rows);
    int i;
    for(i = 0; i < p_cols * p_rows; i++){
        if(p_vals) (*m)->vals[i] = p_vals[i];
        else (*m)->vals[i] = 0.0;
    }
}

void kalman_matrix_init_identity(struct kalman_matrix** m, int p_size){
    *m = (struct kalman_matrix*)malloc(sizeof(struct kalman_matrix));
    (*m)->cols = p_size;
    (*m)->rows = p_size;
    (*m)->vals = (float*)malloc(sizeof(float) * p_size * p_size);
    int i;
    for(i = 0; i < p_size * p_size; i++){
        if(i / p_size == i % p_size) (*m)->vals[i] = 1.0;
        else (*m)->vals[i] = 0.0;
    }
}

void kalman_matrix_init_zero(struct kalman_matrix** m, int p_cols, int p_rows){
    *m = (struct kalman_matrix*)malloc(sizeof(struct kalman_matrix));
    (*m)->cols = p_cols;
    (*m)->rows = p_rows;
    (*m)->vals = (float*)malloc(sizeof(float) * p_cols * p_rows);
    int i;
    for(i = 0; i < p_cols * p_rows; i++){
        (*m)->vals[i] = 0.0;
    }
}

void kalman_matrix_destruct(struct kalman_matrix** m){
    free((*m)->vals);
    free(*m);
    *m = NULL;
}

int kalman_matrix_add(struct kalman_matrix* res,
                      struct kalman_matrix* m1,
                      struct kalman_matrix* m2){
    int i;
    
#ifndef KALMAN_MATRIX_FEELING_LUCKY
    if(!(m1->cols == m2->cols && m1->rows == m2->rows &&
         res->cols == m1->cols && res->rows == m1->rows)){
        return KALMAN_MATRIX_SIZE_MISMATCH;
    }
#endif

    for(i = 0; i < m1->cols * m1->rows; i++){
        res->vals[i] = m1->vals[i] + m2->vals[i];        
    }

    return KALMAN_MATRIX_OK;
}

int kalman_matrix_sub(struct kalman_matrix* res,
                      struct kalman_matrix* m1,
                      struct kalman_matrix* m2){
    int i;

#ifndef KALMAN_MATRIX_FEELING_LUCKY
    if(!(m1->cols == m2->cols && m1->rows == m2->rows &&
         res->cols == m1->cols && res->rows == m1->rows)){
        return KALMAN_MATRIX_SIZE_MISMATCH;
    }
#endif

    for(i = 0; i < m1->cols * m1->rows; i++){
        res->vals[i] = m1->vals[i] - m2->vals[i];        
    }
    return KALMAN_MATRIX_OK;
}

int kalman_matrix_multi(struct kalman_matrix* res,
                        struct kalman_matrix* m1,
                        struct kalman_matrix* m2){
    int i, j, k;

#ifndef KALMAN_MATRIX_FEELING_LUCKY
    if(!(m1->cols == m2->rows &&
         res->cols == m2->cols && res->rows == m1->rows)){
        return KALMAN_MATRIX_SIZE_MISMATCH;
    }
#endif

    for(i = 0; i < m2->cols; i++){
        for(j = 0; j < m1->rows; j++){
            res->vals[i+j*m2->cols] = 0.0;
            for(k = 0; k < m1->cols; k++){
                res->vals[i+j*m2->cols] +=
                    m1->vals[k+j*m1->cols] * m2->vals[i+k*m2->cols];
            }
        }
    }
    return KALMAN_MATRIX_OK;
}

int kalman_matrix_multi_trans(struct kalman_matrix* res,
                              struct kalman_matrix* m1,
                              struct kalman_matrix* m2){
    int i, j, k;

#ifndef KALMAN_MATRIX_FEELING_LUCKY
    if(!(m1->cols == m2->cols &&
         res->cols == m2->rows && res->rows == m1->rows)){
        return KALMAN_MATRIX_SIZE_MISMATCH;
    }
#endif

    for(i = 0; i < m2->rows; i++){
        for(j = 0; j < m1->rows; j++){
            res->vals[i+j*m2->rows] = 0.0;
            for(k = 0; k < m1->cols; k++){
                res->vals[i+j*m2->rows] +=
                    m1->vals[k+j*m1->cols] * m2->vals[k+i*m2->cols];
            }
        }
    }
    return KALMAN_MATRIX_OK;
}

int kalman_matrix_scalar_multi(struct kalman_matrix* res,
                               struct kalman_matrix* m,
                               float p_factor){
    int i;

#ifndef KALMAN_MATRIX_FEELING_LUCKY
    if(!(res->cols == m->cols && res->rows == m->rows)){
        return KALMAN_MATRIX_SIZE_MISMATCH;
    }
#endif

    for(i = 0; i < m->cols * m->rows; i++){
        res->vals[i] = m->vals[i] * p_factor;
    }
    return KALMAN_MATRIX_OK;
}

int kalman_matrix_diag_add(struct kalman_matrix* m, struct kalman_matrix* vec){
    int i;

#ifndef KALMAN_MATRIX_FEELING_LUCKY
    if(!(m->cols == m->rows) || !(vec->cols == 1 && vec->rows == m->rows)){
        return KALMAN_MATRIX_SIZE_MISMATCH;
    }
#endif

    for(i = 0; i < m->cols; i++){
        m->vals[i*(m->cols+1)] += vec->vals[i];
    }
    return KALMAN_MATRIX_OK;
}

int kalman_matrix_trans(struct kalman_matrix* res, struct kalman_matrix* m){

    int i, j;

#ifndef KALMAN_MATRIX_FEELING_LUCKY
    if(!(res->cols == m->cols && res->rows == m->rows)){
        return KALMAN_MATRIX_SIZE_MISMATCH;
    }
#endif

    for(i = 0; i < m->cols; i++){
        for(j = 0; j < m->rows; j++){
            *val(res, j, i) = *val(m, i, j);
        }
    }
    return KALMAN_MATRIX_OK;
}

int kalman_matrix_invert_3x3(struct kalman_matrix* res,
                             struct kalman_matrix* m){
#ifndef KALMAN_MATRIX_FEELING_LUCKY
    if(!(m->cols == 3 && m->rows == 3 && res->cols == 3 && res->rows == 3)){
        return KALMAN_MATRIX_SIZE_MISMATCH;
    }
#endif

    float det = 0.0;

    res->vals[0] = m->vals[4] * m->vals[8] - m->vals[5] * m->vals[7];
    res->vals[3] = m->vals[5] * m->vals[6] - m->vals[3] * m->vals[8];
    res->vals[6] = m->vals[3] * m->vals[7] - m->vals[4] * m->vals[6];

    det = m->vals[0] * res->vals[0] +
          m->vals[1] * res->vals[3] +
          m->vals[2] * res->vals[6];
  
    if(det == 0.0f){
        return KALMAN_MATRIX_NOT_INVERTIBLE;
    }

    res->vals[0] /= det;
    res->vals[3] /= det;
    res->vals[6] /= det;
    res->vals[1] = (m->vals[2] * m->vals[7] - m->vals[1] * m->vals[8]) / det;
    res->vals[2] = (m->vals[1] * m->vals[5] - m->vals[2] * m->vals[4]) / det;
    res->vals[4] = (m->vals[0] * m->vals[8] - m->vals[2] * m->vals[6]) / det;
    res->vals[5] = (m->vals[2] * m->vals[3] - m->vals[0] * m->vals[5]) / det;
    res->vals[7] = (m->vals[1] * m->vals[6] - m->vals[0] * m->vals[7]) / det;
    res->vals[8] = (m->vals[0] * m->vals[4] - m->vals[1] * m->vals[3]) / det;

    return KALMAN_MATRIX_OK;
}

float* val(struct kalman_matrix* m, int p_col, int p_row){
    return m->vals + (p_col + p_row * m->cols);
}

int kalman_matrix_get_cols(struct kalman_matrix* m){
    return m->cols;
}

int kalman_matrix_get_rows(struct kalman_matrix* m){
    return m->rows;
}

float* kalman_matrix_get_vals(struct kalman_matrix* m){
    return m->vals;
}

void kalman_matrix_debug_print(struct kalman_matrix* m){
    int i, j;

    for(j = 0; j < m->rows; j++){
        for(i = 0; i < m->cols; i++){
            printf("%7.3f ", m->vals[i+j*m->cols]);
        }
        printf("\n");
    }
} 

