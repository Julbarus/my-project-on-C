#include "s21_string.h"

s21_size_t s21_strspn(const char *str1, const char *str2) {
    s21_size_t len = 0;
    int lenStr2 = s21_strlen(str2);
    if ((str1 != s21_NULL) && (str2 != s21_NULL)) {  // Валидация
        for (; str1[len] != '\0'; len++) {
            int f = 0;
            for (int j = 0; j < lenStr2; j++) {
                if (str1[len] != str2[j]) {
                    f++;
                } else {
                    break;
                }
            }
            if (f == lenStr2) break;
        }
    }
    return len;
}
