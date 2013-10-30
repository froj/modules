#ifndef KALMAN_MATRIX_H
#define KALMAN_MATRIX_H

struct kalman_matrix;

void kalman_matrix_init(struct kalman_matrix** m,
                        int p_cols,
                        int p_rows,
                        float* p_vals);

void kalman_matrix_init_identity(struct kalman_matrix** m, int p_size);
void kalman_matrix_init_zero(struct kalman_matrix** m, int p_cols, int p_rows);

void kalman_matrix_destruct(struct kalman_matrix** m);

int kalman_matrix_add(struct kalman_matrix* res,
                       struct kalman_matrix* m1,
                       struct kalman_matrix* m2);

int kalman_matrix_sub(struct kalman_matrix* res,
                       struct kalman_matrix* op,
                       struct kalman_matrix* sub);

int kalman_matrix_multi(struct kalman_matrix* res,
                         struct kalman_matrix* m1,
                         struct kalman_matrix* m2);

int kalman_matrix_multi_trans(struct kalman_matrix* res,
                               struct kalman_matrix* m1,
                               struct kalman_matrix* m2);

int kalman_matrix_scalar_multi(struct kalman_matrix* res,
                                struct kalman_matrix* m,
                                float p_factor);

int kalman_matrix_trans(struct kalman_matrix* res,
                         struct kalman_matrix* m);

/** this function will change the matrix! */
int kalman_matrix_diag_add(struct kalman_matrix* m,
                            struct kalman_matrix* vec);

int kalman_matrix_invert_3x3(struct kalman_matrix* res,
                             struct kalman_matrix* m);

//float kalman_matrix_det(struct kalman_matrix* m);

int kalman_matrix_get_cols(struct kalman_matrix* m);
int kalman_matrix_get_rows(struct kalman_matrix* m);
float* kalman_matrix_get_vals(struct kalman_matrix* m);

void kalman_matrix_debug_print(struct kalman_matrix* m);

#define     KALMAN_MATRIX_OK              0
#define     KALMAN_MATRIX_SIZE_MISMATCH   2
#define     KALMAN_MATRIX_NOT_INVERTIBLE  3
#define     KALMAN_MATRIX_NOK             4

#endif
