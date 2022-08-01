#include "s21_string.h"

void* s21_memmove(void *dest, const void *src, s21_size_t n) {
    char *TmpDest = (char *)dest;
    const char *TmpSrc =(const char *)src;
    if ((TmpDest != s21_NULL) && (TmpSrc != s21_NULL)) {  // валидность
        while (n--) {
            *TmpDest++ = *TmpSrc++;
        }
    }
    return dest;
}
