#ifndef SRC_S21_STRING_H_
#define SRC_S21_STRING_H_

#include <stdlib.h>
#include <stdarg.h>
#include <limits.h>
#include <math.h>

typedef long unsigned s21_size_t;

#ifndef s21_NULL
#define s21_NULL ((void*)0)
#endif

#define MAXLEN 150

struct errlist {
    char* name_err;
    int err_num;
    char* describe_err;
};

void  processChar(char* buffer, int* index, const int flags[], va_list va);
void processStr(char* buffer, int* index, int flags[], va_list va);

long double acc_up(long double fraction, int accuracy);
char outputPlusOrMinus(long double *value, int* count);
long double processG(long double value, int flags[], int* count, char* plusOrMinus);
void outputFloat(char *buffer, long double fraction, int *index, int flags[]);
void expExp(char* tempBuff, int* tempIndex, const int flags[], char plusOrMinus, int count);
void zapolnenieFloatInt(char* tempBuff, int* tempIndex, int flags[], long double value);
void processFloatAndExp(char* buffer, int* index, int flags[], va_list va);

void processInt(char* buffer, int* index, int flags[], va_list va);
void zapolneniTempBufferaInt(char* buffer, int* index, int flags[], long long int value);
void processUint(char* buffer, int* index, int flags[], va_list va);
void zapolneniTempBufferaU(char* buffer, int* index, int flags[], long long unsigned int value);
void zapolnenie(char* buffer, int* index, int flags[], const char* tempBuff, const int *tempIndex);

void processPointer(char* buffer, int* index, int flags[], void *value);
void processOctal(char* buffer, int* index, int flags[], va_list va);
void processNum_16(char* buffer, int* index, int flags[], va_list va);
void processN_Spec(int *value, const int *index);

int amountNumber(unsigned long long int n);
void otstupSpase(char* buffer, int* index, int indexSpace, char c);
char* Accuracy_Place(char* ptr_format, int flags[], va_list va);
char* findFlags(int flags[], int* indexForBuffer, char* ptr_format, va_list va);
void specificator(char* buffer, int* indexForBuffer, int flags[], va_list va);
int parser(char* buffer, char* format, va_list va);
int s21_sprintf(char* buffer, char* format, ...);

void at_exit(char *masssage);
long double s21_powl(long double a, int n);

const char* processCharScn(int* amount, const char* ptr_str, const int flags[], va_list va, int* flagCount);
const char* processStrScn(int* amount, const char* ptr_str, int flags[], va_list va, int* flagCount);
void getLenghScn(int flags[], va_list va, unsigned long value);
int convertScn(int value);
int validScn(int num, const char* ptr_str);
const char* getValueNum_8_or_16Scn(int flags[], const char* ptr_str,
    unsigned long* value, int* amount, int num, int* flagCount);
const char* processNum_8_or_16Scn(int* amount, const char* ptr_str, int flags[], va_list va, int* flagCount);
const char* getValueLongLongScn(const char* ptr_str, int flags[],
    long long* value, int* amount, int* flagCount);
const char* processD_and_UScn(int* amount, const char* ptr_str, int flags[], va_list va, int* flagCount);
const char* processIScn(int* amount, const char* ptr_str, int flags[], va_list va, int* flagCount);
const char* processPointerScn(int* amount, const char* ptr_str, int flags[], va_list va, int* flagCount);
const char* getFractionPartScn(const char* ptr_str, int flags[], long double* value);
const char* processFloatAndExpScn(int* amount, const char* ptr_str, int flags[], va_list va, int* flagCount);
const char* specificatorScn(int* amount, const char* ptr_str, int flags[], va_list va, int *flagCount);
const char* widthScn(const char* ptr_format, int flags[]);
const char* findFlagsScn(int flags[], const char* ptr_format);
int parserScn(const char* str, const char* format, va_list va);

int s21_sscanf(const char* str, const char* format, ...);

int def_errlist_mac(struct errlist* errlist_for_mac);
int def_errlist_linux(struct errlist* errlist_for_linux);

void* s21_memchr(const void *str, int c, s21_size_t n);
int s21_memcmp(const void *str1, const void *str2, s21_size_t n);
void* s21_memcpy(void *dest, const void *src, s21_size_t n);
void* s21_memmove(void *dest, const void *src, s21_size_t n);
void *s21_memset(void *str, int c, s21_size_t n);
char *s21_strcat(char *dest, const char *src);
char *s21_strncat(char *dest, const char *src, s21_size_t n);
char *s21_strchr(const char  *str, int c);
int s21_strcmp(const char *str1, const char *str2);
int s21_strncmp(const char *str1, const char *str2, s21_size_t n);
char* s21_strcpy(char *dest, const char *src);
char* s21_strncpy(char *dest, const char *src, s21_size_t n);
s21_size_t s21_strcspn(const char *str1, const char *str2);
char* s21_strerror(int errnum);
s21_size_t s21_strlen(const char *str);
char *s21_strpbrk(const char *str1, const char *str2);
char *s21_strrchr(const char *src, int c);
s21_size_t s21_strspn(const char *str1, const char *str2);
char *s21_strstr(const char* haystack, const char* needle);
char *s21_strtok(char *str, const char *delim);
void *s21_to_upper(const char *str);
void *s21_to_lower(const char *str);
void *s21_insert(const char *src, const char *str, s21_size_t start_index);
void *s21_trim(const char *src, const char *trim_chars);


#endif  // SRC_S21_STRING_H_
