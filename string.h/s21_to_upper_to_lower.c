#include "s21_string.h"

void *s21_to_upper(const char *str) {
    char *buff = s21_NULL;

    if (str) {
        s21_size_t length = s21_strlen(str);
        buff = (char*)calloc((length + 1), sizeof(char));
        for (s21_size_t i = 0; i <= length; i++)
            (str[i] >= 97 && str[i] <= 122) ? (buff[i] = str[i] - 32) : (buff[i] = str[i]);
    }
    return buff;
}

void *s21_to_lower(const char *str) {
    char *buff = s21_NULL;

    if (str) {
        s21_size_t length = s21_strlen(str);
        buff = (char*)malloc(sizeof(char) * (length + 1));
        for (s21_size_t i = 0; i <= length; i++)
                (str[i] >= 65 && str[i] <= 90) ? (buff[i] = str[i] + 32) : (buff[i] = str[i]);
    }
    return buff;
}
