#include "../s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  if (rows <= 0 || columns <= 0 || result == NULL)
    return ERROR_INCORRECT_MATRIX;
  int resulting_code = OK;

  result->matrix = (double **)malloc(rows * sizeof(double *));
  if (result->matrix == NULL) {
    s21_remove_matrix(result);
    resulting_code = ERROR_INCORRECT_MATRIX;
  }
  if (resulting_code == OK) {
    for (int i = 0; i < rows && resulting_code == OK; i++) {
      result->matrix[i] = (double *)malloc(columns * sizeof(double));
      if (result->matrix[i] == NULL) {
        s21_remove_matrix(result);
        resulting_code = ERROR_INCORRECT_MATRIX;
      }
    }
  }
  if (resulting_code == OK) {
    result->rows = rows;
    result->columns = columns;
    null_out_matrix(result);
  }

  return resulting_code;
}

void s21_remove_matrix(matrix_t *A) {
  for (int i = 0; i < A->rows; i++) {
    if (A->matrix[i]) free(A->matrix[i]);
  }
  if (A->matrix) free(A->matrix);
  A->matrix = NULL;
}

/**
 * @brief Checks if matrices are equal
 *
 * @note
 *
 * @param A First matrix
 * @param B Second matrix
 *
 * @return SUCCESS (1) if matrices are equal, FAILURE (0) otherwise
 */
int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  if (A == NULL || B == NULL) return ERROR_INCORRECT_MATRIX;
  if (NOT(matrices_are_the_same_size(A, B))) return FAILURE;

  int resulting_code = SUCCESS;
  for (int i = 0; i < A->rows && resulting_code == SUCCESS; i++) {
    for (int j = 0; j < A->columns && resulting_code == SUCCESS; j++) {
      if (fabs(A->matrix[i][j] - B->matrix[i][j]) >= _maximum_fault)
        resulting_code = FAILURE;
    }
  }
  return resulting_code;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A == NULL || B == NULL || result == NULL) return ERROR_INCORRECT_MATRIX;
  if (NOT(matrices_are_the_same_size(A, B))) return CALCULATION_ERROR;

  int resulting_code = OK;
  null_out_matrix(result);
  for (int i = 0; i < A->rows && resulting_code == OK; i++) {
    for (int j = 0; j < A->columns && resulting_code == OK; j++) {
      double sum = A->matrix[i][j] + B->matrix[i][j];
      if (isinf(sum))
        resulting_code = CALCULATION_ERROR;
      else
        result->matrix[i][j] = sum;
    }
  }
  return resulting_code;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A == NULL || B == NULL || result == NULL) return ERROR_INCORRECT_MATRIX;
  if (NOT(matrices_are_the_same_size(A, B))) return CALCULATION_ERROR;

  int resulting_code = OK;
  null_out_matrix(result);
  for (int i = 0; i < A->rows && resulting_code == OK; i++) {
    for (int j = 0; j < A->columns && resulting_code == OK; j++) {
      double sub = A->matrix[i][j] - B->matrix[i][j];
      if (isinf(sub))
        resulting_code = CALCULATION_ERROR;
      else
        result->matrix[i][j] = sub;
    }
  }
  return resulting_code;
}

/**
 * @brief The product of the matrix A = m × n by the number λ is the matrix B =
 * m × n = λ × A whose elements are defined by the equations B = λ × A(i,j).
 *
 * @param A First matrix
 * @param number Number
 * @param result Resulting matrix
 *
 * @return Matrix
 */
int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  if (A == NULL || result == NULL) return ERROR_INCORRECT_MATRIX;

  int resulting_code = OK;
  null_out_matrix(result);
  for (int i = 0; i < A->rows && resulting_code == OK; i++) {
    for (int j = 0; j < A->columns && resulting_code == OK; j++) {
      double product = A->matrix[i][j] * number;
      if (isinf(product))
        resulting_code = CALCULATION_ERROR;
      else
        result->matrix[i][j] = product;
    }
  }
  return resulting_code;
}

/**
 * @brief The product of A = m × k by B = k × n is a matrix C = m × n = A × B of
 * size m × n whose elements are defined by the equation C(i,j) = A(i,1) ×
 * B(1,j) + A(i,2) × B(2,j) + ... + A(i,k) × B(k,j).
 *
 * @param A First matrix
 * @param B Second matrix
 * @param result Resulting matrix
 *
 * @return Matrix
 */
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A == NULL || B == NULL || result == NULL) return ERROR_INCORRECT_MATRIX;
  if (A->columns != B->rows) return CALCULATION_ERROR;

  int resulting_code = OK;
  null_out_matrix(result);
  for (int i = 0; i < A->rows && resulting_code == OK; i++) {
    for (int j = 0; j < B->columns && resulting_code == OK; j++) {
      result->matrix[i][j] = 0;
      for (int k = 0; k < A->columns && resulting_code == OK; k++) {
        result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        if (isinf(result->matrix[i][j])) resulting_code = CALCULATION_ERROR;
      }
    }
  }
  return resulting_code;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  if (A == NULL || result == NULL) return ERROR_INCORRECT_MATRIX;

  int resulting_code = OK;
  null_out_matrix(result);
  for (int i = 0; i < A->rows && resulting_code == OK; i++) {
    for (int j = 0; j < A->columns && resulting_code == OK; j++) {
      result->matrix[j][i] = A->matrix[i][j];
    }
  }
  return resulting_code;
}

/**
 * @brief Calculates the determinant of a matrix into the result
 *
 * @note Minor M(i,j) is a (n-1)-order determinant obtained by deleting out the
 * i-th row and the j-th column from the matrix A.\
 *
 * The algebraic complement of a matrix element is the value of the minor
 * determinants multiplied by -1^(i+j).
 *
 * Hierarchy of functions:
 * s21_calc_complements() {
 *   matrix_minor();
 *   s21_determinant() {
 *     triangulation_in_place() {
 *       permutation_of_rows_in_place();
 *     }
 *   }
 * }
 *
 * @param A First matrix
 * @param result Resulting matrix
 *
 * @return Matrix
 */
int s21_calc_complements(matrix_t *A, matrix_t *result) {
  if (A == NULL || result == NULL) return ERROR_INCORRECT_MATRIX;
  if (A->rows != A->columns) return ERROR_INCORRECT_MATRIX;

  int resulting_code = OK;
  double determinant = 0;
  null_out_matrix(result);
  for (int i = 0; i < A->rows && resulting_code == OK; i++) {
    for (int j = 0; j < A->columns && resulting_code == OK; j++) {
      matrix_t minor;
      resulting_code = s21_create_matrix(A->rows, A->columns, &minor);
      resulting_code = matrix_minor(A, &minor, i, j);
      result->matrix[i][j] =
          s21_determinant(&minor, &determinant) * pow(-1, i + j);
      s21_remove_matrix(&minor);
    }
  }
  return resulting_code;
}

int matrix_minor(matrix_t *A, matrix_t *minor, int i, int j) {
  if (A == NULL || minor == NULL) return ERROR_INCORRECT_MATRIX;

  int resulting_code = OK;
  resulting_code = s21_create_matrix(A->rows - 1, A->columns - 1, minor);
  if (resulting_code == OK) {
    for (int k = 0; k < A->rows - 1; k++) {
      for (int l = 0; l < A->columns - 1; l++) {
        minor->matrix[k][l] = A->matrix[k < i ? k : k + 1][l < j ? l : l + 1];
      }
    }
  }
  return resulting_code;
}

/**
 * @brief Calculates the determinant of a matrix into the result
 *
 * @note Finds the determinant by Gauss elimination method.
 *
 * @param A First matrix
 * @param result Resulting matrix
 *
 * @return Matrix
 *
 */
int s21_determinant(matrix_t *A, double *result) {
  if (A == NULL || result == NULL) return ERROR_INCORRECT_MATRIX;
  if (A->rows != A->columns) return CALCULATION_ERROR;

  int resulting_code = OK;
  double determinant = 1;
  matrix_t temp;
  s21_create_matrix(A->rows, A->columns, &temp);
}

/**
 * @brief Triangulates a matrix in place
 *
 * @note To "zero" the elements of the i-th column of the matrix, it is
 * sufficient to add the i-th row multiplied by -a[j][i]/a[i][i] to all rows
 * with numbers j = i+1, ... n.
 *
 * Example:
 * {{10, 2, 3},
 * {4, 5, 6},
 * {7, 8, 9}}
 *
 * 0th column: j = 1; i = 0 => {4 - 10 * 4/10, 5 - 10 * 4/10, 6 - 10 * 4/10} =
 * {0, 1, 2}
 * j = 2; i = 0 => {7 - 10 * 7/10, 8 - 10 * 7/10, 9 - 10 * 7/10} =
 * {0, 1, 2}
 *
 * {{10, 2, 3},
 * {0, 1, 2},
 * {0, 1, 2}}
 *
 * 1st column: j = 2; i = 1 => {0 - 0 * 1/1, 1 - 1 * 1/1, 2 - 2 * 1/1} =
 * {0, 0, 0}
 *
 * {{10, 2, 3},
 * {0, 1, 2},
 * {0, 0, 0}}
 *
 * Then, the determinant equals 10 * 1 * 0 = 0
 *
 * When performing such an operation, division by
 * zero may occur if the element on the main diagonal is equal to zero - in this
 * case, the matrix rows are permuted.
 *
 * The most effective approach to permuting
 * rows is to permutate the i-th row with the row that has the maximum element
 * in the i-th column. It has been proven that when performing such a
 * permutation for each row (and not only when a[i][i] is equal to zero) for a
 * matrix with a non-zero determinant, we will always find a row for
 * replacement.
 *
 * Source: https://pro-prof.com/forums/topic/matrix-triangulation
 */
int triangulation_in_place(matrix_t *A) {
  if (A == NULL) return ERROR_INCORRECT_MATRIX;

  int resulting_code = OK;
  for (int i = 0; i < A->rows && resulting_code == OK; i++) {
    permutation_of_rows_in_place(A, i);
    for (int j = i + 1; j < A->columns && resulting_code == OK; j++) {
      double multiplier = A->matrix[j][i] / A->matrix[i][i];
      if (isinf(multiplier)) {
        permutation_of_rows_in_place(A, i);
      }
      for (int k = i; k < A->columns && resulting_code == OK; k++) {
        A->matrix[j][k] -= multiplier * A->matrix[i][k];
      }
    }
  }
  return resulting_code;
}

/**
 * @brief Permutates the i-th row with the row that has the maximum element
 * in the i-th column.
 */
void permutation_of_rows_in_place(matrix_t *A, int i) {
  int max_index = i;
  double max_element = A->matrix[i][i];
  for (int j = i + 1; j < A->columns; j++) {
    if (A->matrix[j][i] > max_element) {
      max_element = A->matrix[j][i];
      max_index = j;
    }
  }

  if (max_index != i) {
    double *temp = A->matrix[i];
    A->matrix[i] = A->matrix[max_index];
    A->matrix[max_index] = temp;
  }
}

void null_out_matrix(matrix_t *A) {
  if (A == NULL) return;
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      A->matrix[i][j] = 0;
    }
  }
}

bool matrices_are_the_same_size(matrix_t *A, matrix_t *B) {
  return A->rows == B->rows && A->columns == B->columns;
}