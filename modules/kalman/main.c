#include <stdio.h>
#include "kalman_matrix.h"


int main(void){
    struct kalman_matrix* mat0;
    struct kalman_matrix* mat1;
    struct kalman_matrix* mat2;
    struct kalman_matrix* mat3;
    struct kalman_matrix* mat4;
    
    struct kalman_matrix* buf;
    struct kalman_matrix* res;

    float vals0[] = {1.0, 2.0, 3.0,
                    4.0, 5.0, 6.0,
                    7.0, 8.0, 0.0};

    float vals3[] = {1.0, 0.0, 0.0,
                     0.0, 1.0, 0.0,
                     0.0, 0.0, 1.0,
                     0.0, 0.0, 1.0,
                     0.0, 1.0, 0.0,
                     1.0, 0.0, 0.0};

    float vals4[] = {1.0,
                     2.0,
                     3.0};

    printf("Init matrix (M0):\n");
    kalman_matrix_init(&mat0, 3, 3, vals0);
    kalman_matrix_debug_print(mat0);
    printf("\n");
    printf("Init identity (M1):\n");
    kalman_matrix_init_identity(&mat1, 3);
    kalman_matrix_debug_print(mat1);
    printf("\n");
    printf("Init all zeros (M2):\n");
    kalman_matrix_init_zero(&mat2, 3, 6);
    kalman_matrix_debug_print(mat2);
    printf("\n");
    printf("Init matrix (M3):\n");
    kalman_matrix_init(&mat3, 3, 6, vals3);
    kalman_matrix_debug_print(mat3);
    printf("\n");
    printf("Init matrix (M4):\n");
    kalman_matrix_init(&mat4, 1, 3, vals4);
    kalman_matrix_debug_print(mat4);
    printf("\n");
    printf("BUF\n");
    kalman_matrix_init_zero(&buf, 6, 6);
    kalman_matrix_debug_print(buf);
    printf("\n");
    printf("RES\n");
    kalman_matrix_init_zero(&res, 3, 3);
    kalman_matrix_debug_print(res);
    printf("\n");

    printf("RES = M0 + M1\n");
    kalman_matrix_add(res, mat0, mat1);
    kalman_matrix_debug_print(res);
    printf("\n");

    printf("RES = M0 - M1\n");
    kalman_matrix_sub(res, mat0, mat1);
    kalman_matrix_debug_print(res);
    printf("\n");

    printf("M2 = M3 * M0\n");
    kalman_matrix_multi(mat2, mat3, mat0);
    kalman_matrix_debug_print(mat2);
    printf("\n");

    printf("BUF = M2 * transp(M3)\n");
    kalman_matrix_multi_trans(buf, mat2, mat3);
    kalman_matrix_debug_print(buf);
    printf("\n");

    printf("BUF = BUF * 100\n");
    kalman_matrix_scalar_multi(buf, buf, 100.0f);
    kalman_matrix_debug_print(buf);
    printf("\n");

    printf("RES = transp(M0)\n");
    kalman_matrix_trans(res, mat0);
    kalman_matrix_debug_print(res);
    printf("\n");

    printf("M1 = M1 diagonal_add M4\n");
    kalman_matrix_diag_add(mat1, mat4);
    kalman_matrix_debug_print(mat1);
    printf("\n");

    printf("RES = inv(M0)\n");
    kalman_matrix_invert_3x3(res, mat0);
    kalman_matrix_debug_print(res);
    printf("\n");

//void kalman_matrix_destruct(struct kalman_matrix** m);
//
//
//
///** this function will change the matrix! */
//int kalman_matrix_diag_add(struct kalman_matrix* m,
//                            struct kalman_matrix* vec);
//
//int kalman_matrix_invert_3x3(struct kalman_matrix* res,
//                             struct kalman_matrix* m);
//
}
