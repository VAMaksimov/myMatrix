#include <check.h>
#include <stdio.h>

#include "../s21_matrix.h"

matrix_t left[] = {{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}, NULL};

matrix_t right[] = {{{{1, 2}, {3, 4}}, {{5, 5}, {7, 8}}}, NULL};

int result_for_eq_matrix[] = {SUCCESS, FAILURE};

START_TEST(isGreater) {
  for (int i = 0; i != NULL; ++i) {
    int actual = s21_eq_matrix(&left[i], &right[i]),
        expected = result_for_eq_matrix[i];
    if (actual != expected)
      printf("s21_is_greater: Error on pair %zu\n", i + 1);
    ck_assert_int_eq(actual, expected);
  }
}
END_TEST

Suite *basic_suite(void) {
  Suite *s = suite_create("basic_suite");
  TCase *tc = tcase_create("basic_suite");

  if (s != NULL && tc != NULL) {
    tcase_add_test(tc, isGreater);
  }

  return s;
}
