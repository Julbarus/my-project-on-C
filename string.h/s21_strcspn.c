#include "s21_string.h"

s21_size_t s21_strcspn(const char *str1, const char *str2) {
    s21_size_t len = 0;
    if ((str1 != s21_NULL) && (str2 != s21_NULL)) {  // Валидация
        int fl = 0;
        for (; str1[len] != '\0'; len++) {
            for (int i = 0; str2[i] != '\0'; i++) {
                if (str1[len] == str2[i]) {
                    fl = 1;
                    break;
                }
            }
            if (fl) break;
        }
    }
    return len;
}
