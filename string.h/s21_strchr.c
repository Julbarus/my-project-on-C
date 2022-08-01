#include "s21_string.h"

char *s21_strchr(const char *str, int c) {
    char *ch = 0;
    int len = s21_strlen(str) + 1;
    for (int i = 0; i < len; i++) {
        if (str[i] == (char)c) {
            ch = (char*)&str[i];
            break;
        }
    }
    return ch;
}

char *s21_strrchr(const char *src, int c) {
    char *r = 0;
    int end = s21_strlen(src) + 1;
    while (0 < end--) {
        if (src[end] == c) {
            r = (char*)&src[end];
            break;
        }
    }
    return r;
}
