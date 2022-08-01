#include "s21_string.h"

char *s21_strpbrk(const char *str1, const char *str2) {
    int i, j, flag = 0;
    for (i = 0; str1[i] != '\0'; i++) {
        for (j = 0; str2[j] != '\0'; j++) {
            if (str1[i] == str2[j])
                flag = 1;
            if (flag == 1)
                break;
        }
        if (flag == 1)
            break;
    }
    return flag == 1 ? ((char *) str1 + i) : 0;
}
