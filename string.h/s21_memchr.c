#include "s21_string.h"

void* s21_memchr(const void *str, int c, s21_size_t n) {
    void *a = s21_NULL;
    if (n != 0 && str != s21_NULL) {
        const unsigned char *p = str;
        for (s21_size_t i = 0; i < n; i++) {
            if ((unsigned char)c == *p++) {
                a = (void*)(p-1);
                break;
            }
        }
    }
    return a;
}
