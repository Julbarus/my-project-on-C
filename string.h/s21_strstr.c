#include "s21_string.h"

char *s21_strstr(const char *str1, const char *str2) {
  char *str1_tmp = (char *)str1;
  char *str2_tmp = (char *)str2;
  char *res = s21_NULL;
  if (*str2_tmp == '\0') res = str1_tmp;
  for (s21_size_t i = 0; i < s21_strlen(str1_tmp); i++) {
    if (str1_tmp[i] == *str2_tmp) {
      char *tmp = &str1_tmp[i];
      int tmp_res = s21_strncmp(tmp, str2_tmp, s21_strlen(str2_tmp));
      if (tmp_res == 0) {
        res = &str1_tmp[i];
        break;
      }
    }
  }
  return res;
}
