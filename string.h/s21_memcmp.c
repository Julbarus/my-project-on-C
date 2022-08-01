#include "s21_string.h"

int s21_memcmp(const void *str1, const void *str2, s21_size_t n) {
    const char *s1 = str1;
    const char *s2 = str2;
    int flag = 0;
    if (n != 0) {
            do {
                if (*s1 != *s2) {
                    flag = 1;
                    break;
                }
                s1++;
                s2++;
                n--;
            } while (n != 0);
    }
    return (flag == 0 ? 0: *s1 - *s2);
}

