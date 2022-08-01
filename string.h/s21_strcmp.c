#include "s21_string.h"

int s21_strcmp(const char *str1, const char *str2) {
    int answer = 1;
    if (str1 != s21_NULL && str2 != s21_NULL) {
        while (*str1 == *str2) {
            if (*str1++ == '\0') {
                answer =  0;
                break;
            }
            str2++;
        }
    } else {
        answer = 0;
    }
    return (answer == 0 ? 0 : (int)(*str1 - *str2));
}

int s21_strncmp(const char *str1, const char *str2, s21_size_t n) {
    int answer = 1;
    if (n == 0 || str1 == s21_NULL || str2 == s21_NULL) {
        answer = 0;
    } else {
        while (*str1 == *str2 && 0 < --n) {
            if (*str1++ == '\0') {
                answer = 0;
                break;
            }
            str2++;
        }
    }
    return (answer == 0 ? 0 : (int)(*str1 - *str2));
}
