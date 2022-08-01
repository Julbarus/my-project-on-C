#include "s21_string.h"

char* s21_strcpy(char *dest, const char *src) {
    int i;
    for (i = 0; src[i] != '\0'; i++) {
        dest[i] = src[i];
    }
    dest[i] = '\0';
    return dest;
}

char* s21_strncpy(char *dest, const char *src, s21_size_t n) {
    s21_size_t i = 0;
    while (i < n) {
        dest[i] = src[i];
        i++;
    }

    return dest;
}
