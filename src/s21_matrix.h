#ifndef S21_DECIMAL_H
#define S21_DECIMAL_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define bool _Bool
#define true 1
#define false 0

#define CONVERTATION_SUCCESS 0
#define CONVERTATION_ERROR 1

#define ADD_OK 0
#define NUMBER_TOO_LARGE 1
#define NUMBER_TOO_SMALL 2
#define DIVISION_BY_ZERO 3

#ifndef NULL
#define NULL ((void *)0)
#endif

typedef struct {
  int bits[7];
} s21_long_decimal;

typedef struct {
  int bits[4];
} s21_decimal;

typedef union {
  int ui;
  float fl;
} fbits;

#define ROW_NUMBER (int)(sizeof(s21_decimal) / sizeof(int))
#define COLUMN_NUMBER (int)(sizeof(int) * 8)
#define DECIMAL_SIZE (ROW_NUMBER * COLUMN_NUMBER)
#define ADDITIONAL_INT_BIT (COLUMN_NUMBER - 1)
#define VALUE_PART_SIZE (3 * COLUMN_NUMBER)  // 96
#define MINUS_BIT_INDEX (ROW_NUMBER * COLUMN_NUMBER - 1)

#define LONG_ROW_NUMBER (int)(sizeof(s21_long_decimal) / sizeof(int))
#define LONG_VALUE_PART_SIZE (6 * COLUMN_NUMBER)
#define LONG_MINUS_BIT_INDEX (LONG_ROW_NUMBER * COLUMN_NUMBER - 1)

#define __INT_MIN__ 0x80000000

#define max(x, y) ((x) > (y) ? (x) : (y))
#define abs(x) ((x) < 0 ? -(x) : (x))
#define getRow(x) (x / COLUMN_NUMBER)
#define getColumn(x) (x % COLUMN_NUMBER)

// minor_utils
void resetDecimal(s21_decimal *dst);
void nullOutDecimal(s21_decimal *dst);
void nullOutLongDecimal(s21_long_decimal *dst);

void negateDecimal(s21_decimal *dst);

bool isSetBit(s21_decimal number, int index);
bool isSetLongBit(s21_long_decimal dst, int index);
bool isScaleSet(s21_decimal dst);
void setBit(s21_decimal *dst, int index);
void setLongBit(s21_long_decimal *dst, int index);

void inverseBit(s21_decimal *dst, int index);

void resetBit(s21_decimal *dst, int index);
void resetLongBit(s21_long_decimal *dst, int index);
bool isDecimalZero(s21_decimal num);

int getScale(s21_decimal num);
void setScale(s21_decimal *num, int scale);
bool isCorrectDecimal(s21_decimal *num);

void alignScale(s21_decimal *value_1, s21_decimal *value_2, int *errorType);
void multiplyBy10(s21_decimal *src, int *errorType);
void shift_left(s21_decimal *src, int shift, int *errorType);

bool isIntPartBiggerOrEqual(s21_decimal value_1, s21_decimal value_2);
void copySign(s21_decimal value_1, s21_decimal *result);
int determineTheSizeDifference(s21_decimal value_1, s21_decimal value_2);

void printDecimal(s21_decimal value);
int div_10(s21_decimal *value);

// convertation
// int countSignificantDigits(float src);
// int calculate_scale(float src);
int getFloatExp(const float *value);

s21_decimal *setBitFloat(s21_decimal *value, int pos, int bit);
void setNegativeSign(s21_decimal *value, int bit);

int s21_from_int_to_decimal(int src, s21_decimal *dst);
int s21_from_decimal_to_int(s21_decimal src, int *dst);
int s21_from_decimal_to_float(s21_decimal src, float *dst);
int s21_from_float_to_decimal(float src, s21_decimal *dst);

// arithmetics
int s21_mul(s21_decimal value_1, s21_decimal value_2, s21_decimal *result);
int s21_add(s21_decimal value_1, s21_decimal value_2, s21_decimal *result);
int s21_sub(s21_decimal value_1, s21_decimal value_2, s21_decimal *result);
int s21_div(s21_decimal value_1, s21_decimal value_2, s21_decimal *result);

int addWrapper(s21_decimal value_1, s21_decimal value_2, s21_decimal *result);
void addition(s21_decimal a, s21_decimal b, s21_decimal *result,
              int *errorType);
void subtraction(s21_decimal a, s21_decimal b, s21_decimal *result,
                 int *errorType);
// comparison
int s21_is_greater(s21_decimal a, s21_decimal b);
int s21_is_less(s21_decimal a, s21_decimal b);
int s21_is_equal(s21_decimal a, s21_decimal b);
int s21_is_not_equal(s21_decimal a, s21_decimal b);
int s21_is_greater_or_equal(s21_decimal a, s21_decimal b);
int s21_is_less_or_equal(s21_decimal a, s21_decimal b);

// others
int s21_truncate(s21_decimal value, s21_decimal *result);
int s21_floor(s21_decimal value, s21_decimal *result);
int s21_round(s21_decimal value, s21_decimal *result);
int s21_negate(s21_decimal value, s21_decimal *result);
int div_10(s21_decimal *value);

#endif  // S21_DECIMAL_H
