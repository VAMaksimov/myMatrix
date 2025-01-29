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
      if (A->matrix[i][j] != B->matrix[i][j]) resulting_code = FAILURE;
    }
  }
  return resulting_code;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A == NULL || B == NULL || result == NULL) return ERROR_INCORRECT_MATRIX;
  if (NOT(matrices_are_the_same_size(A, B))) return CALCULATION_ERROR;

  int resulting_code = OK;
  for (int i = 0; i < A->rows && resulting_code == OK; i++) {
    for (int j = 0; j < A->columns && resulting_code == OK; j++) {
      result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    }
  }
  return resulting_code;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A == NULL || B == NULL || result == NULL) return ERROR_INCORRECT_MATRIX;
  if (NOT(matrices_are_the_same_size(A, B))) return CALCULATION_ERROR;

  int resulting_code = OK;
  for (int i = 0; i < A->rows && resulting_code == OK; i++) {
    for (int j = 0; j < A->columns && resulting_code == OK; j++) {
      result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
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
  for (int i = 0; i < A->rows && resulting_code == OK; i++) {
    for (int j = 0; j < A->columns && resulting_code == OK; j++) {
      result->matrix[i][j] = A->matrix[i][j] * number;
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
  for (int i = 0; i < A->rows && resulting_code == OK; i++) {
    for (int j = 0; j < B->columns && resulting_code == OK; j++) {
      result->matrix[i][j] = 0;
      for (int k = 0; k < A->columns && resulting_code == OK; k++) {
        result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
      }
    }
  }
  return resulting_code;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  if (A == NULL || result == NULL) return ERROR_INCORRECT_MATRIX;

  int resulting_code = OK;
  for (int i = 0; i < A->rows && resulting_code == OK; i++) {
    for (int j = 0; j < A->columns && resulting_code == OK; j++) {
      result->matrix[j][i] = A->matrix[i][j];
    }
  }
  return resulting_code;
}

/**
 * @brief Calculates the determinant of a matrix
 *
 * @note Minor M(i,j) is a (n-1)-order determinant obtained by deleting out the
 * i-th row and the j-th column from the matrix A.\
 *
 * The algebraic complement of a matrix element is the value of the minor
 * determinants multiplied by -1^(i+j).
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
  for (int i = 0; i < A->rows && resulting_code == OK; i++) {
    for (int j = 0; j < A->columns && resulting_code == OK; j++) {
      matrix_t minor;
      resulting_code = s21_create_matrix(A->rows, A->columns, &minor);
      minor = matrix_minor(A, i, j);
      result->matrix[i][j] =
          s21_determinant(&minor, &determinant) * pow(-1, i + j);
      s21_remove_matrix(&minor);
    }
  }
  return resulting_code;
}

matrix_t matrix_minor(matrix_t *A, int i, int j) {
  matrix_t minor;
  s21_create_matrix(A->rows - 1, A->columns - 1, &minor);
  for (int k = 0; k < A->rows - 1; k++) {
    for (int l = 0; l < A->columns - 1; l++) {
      minor.matrix[k][l] = A->matrix[k < i ? k : k + 1][l < j ? l : l + 1];
    }
  }
  return minor;
}

/**
 * @brief Calculates the determinant of a matrix
 *
 * @note Finds the determinant by Gauss elimination method.
 * At the first stage, the so-called direct move is carried out, when by means
 * of elementary transformations over the rows the system is reduced to a
 * stepped or triangular form, or it is established that the system is
 * incompatible. To do this, a non-zero is selected from the elements of the
 * first column of the matrix, the row containing it is moved to the uppermost
 * position, making this row the first. Then the non-zero elements of the first
 * column of all the underlying rows are zeroed by subtracting from each row the
 * first row multiplied by the ratio of the first element of these rows to the
 * first element of the first row. After the specified transformations have been
 * performed, the first row and the first column are mentally crossed out and
 * continued until a matrix of zero size remains. If at any of the iterations a
 * non-zero is not found among the elements of the first column, then they move
 * on to the next column and perform a similar operation. At the second stage,
 * the so-called reverse move is carried out, the essence of which is to express
 * all the resulting basic variables through non-basic ones and construct a
 * fundamental system of solutions, or, if all the variables are basic, then to
 * express in numerical form the only solution of the system of linear
 * equations. This procedure begins with the last equation, from which the
 * corresponding basic variable is expressed (and there is only one) and
 * substituted into the previous equations, and so on, moving up the "steps".
 * Each line corresponds to exactly one basic variable, therefore, at each step,
 * except for the last (the topmost), the situation exactly repeats the case of
 * the last line.
 *
 * source:
 * https://ru.wikipedia.org/wiki/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%93%D0%B0%D1%83%D1%81%D1%81%D0%B0
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
  double determinant = 0;
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