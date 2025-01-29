#ifndef S21_MATRIX_H
#define S21_MATRIX_H

#include <math.h>

#define bool _Bool
#define true 1
#define false 0
#define NOT(x) !x

#ifndef NULL
#define NULL ((void *)0)
#endif

#define _maximum_fault 1e-7
#define fabs(x) (x < 0 ? -x : x)

#define OK 0
#define ERROR_INCORRECT_MATRIX 1
#define CALCULATION_ERROR 2  // mismatched matrix sizes;
// matrix for which calculations cannot be performed, etc.

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);

#define SUCCESS 1
#define FAILURE 0

int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

#endif  // S21_MATRIX_H
