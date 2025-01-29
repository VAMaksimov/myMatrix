#ifndef TEST_H
#define TEST_H

#include <check.h>
#include <stdio.h>

#include "../s21_matrix.h"

void assertMatrix(matrix_t expected, matrix_t actual);
Suite *basic_suite(void);

#endif  // TEST_H