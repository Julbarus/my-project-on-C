#include "s21_string.h"

char *s21_strcat(char *dest, const char *src) {
    char *s = dest;
    while (*dest) {
        dest++;
    }
    while (*src) {
        *dest++ = *src++;
    }
    *dest = '\0';
    return s;
}

char *s21_strncat(char *dest, const char *src, s21_size_t n) {
    char *c = dest;
    while (*dest) {
        dest++;
    }
    while (*src && 0 < n--) {
        *dest++ = *src++;
    }
    *dest = '\0';
    return c;
}
