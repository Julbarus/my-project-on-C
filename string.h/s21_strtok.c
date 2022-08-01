#include "s21_string.h"

char *s21_strtok(char *str, const char *delim) {
    static char *next = s21_NULL;
    int flag_1 = 1, flag_2 = 1;
    if (str) {
      next = str;
    }
    if (next && *next == 0) {
      flag_1 = 0;
    }
    char* c = next;
    if (c && flag_1) {
        while (s21_strchr(delim, *c))
            c++;
    }
    if (c && *c == 0 && flag_1) {
      flag_1 = 0;
    }
    char* word = c;
    if (c && flag_1) {
        while (s21_strchr(delim, *c) == 0)
        c++;
    }
    if (c && *c == 0 && flag_1) {
      next = c;
      flag_2 = 0;
    }
    if (c && flag_1 && flag_2) {
        *c = 0;
        next = c + 1;
    }
    return (flag_1 == 0) ? s21_NULL : word;
}
