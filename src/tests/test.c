#include "test.h"

void assert_matrix(matrix_t expected, matrix_t actual) {
  ck_assert_int_eq(expected.rows, actual.rows);
  ck_assert_int_eq(expected.columns, actual.columns);
  for (int i = 0; i < expected.rows; i++) {
    for (int j = 0; j < expected.columns; j++) {
      if (s21_eq_matrix(&expected, &actual) == FAILURE) {
        printf("ERROR.\nactual:\n");
        print_matrix(actual);
        printf("expected:\n");
        print_matrix(expected);
      }
      ck_assert_double_eq_tol(expected.matrix[i][j], actual.matrix[i][j],
                              _maximum_fault);
    }
  }
}

void print_matrix(matrix_t A) {
  for (int i = 0; i < A.rows; i++) {
    for (int j = 0; j < A.columns; j++) {
      printf("%f ", A.matrix[i][j]);
    }
    printf("\n");
  }
}

int main(void) {
  int fails_count = 0;
  Suite *suite_list[] = {basic_suite(), NULL};

  for (Suite **current_suite = suite_list; *current_suite != NULL;
       current_suite++) {
    SRunner *runner = srunner_create(*current_suite);
    srunner_run_all(runner, CK_NORMAL);
    int suite_failed = srunner_ntests_failed(runner);
    fails_count += suite_failed;
    if (suite_failed)
      printf("\033[31mFAILED\033[0m\n");
    else
      printf("\033[32mPASSED\033[0m\n");
    srunner_free(runner);
  }
  if (fails_count)
    printf("\n\033[31m========= FAILED: %d =========\n", fails_count);
  return 0;
}