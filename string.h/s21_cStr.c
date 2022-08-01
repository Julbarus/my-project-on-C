#include "s21_string.h"

void  processChar(char* buffer, int* index, const int flags[], va_list va) {
    if (flags[0] || flags[1] > 1 || flags[2] || flags[4] != -1) {
        *index = -1;
    } else if ((flags[5] > 0 && flags[5] < 3) || flags[5] == 4) {
        *index = -1;
    } else if (*index != -1) {
        buffer[*index] = (char)va_arg(va, long int);
        *index+=1;
    }
}

void processStr(char* buffer, int* index, int flags[], va_list va) {
    if (flags[0] || flags[1] > 1 || flags[2]) {
        *index = -1;
    } else if ((flags[5] > 0 && flags[5] < 3) || flags[5] == 4) {
        *index = -1;
    } else if (*index != -1) {
        if (flags[5] == 3) {
            wchar_t *value = 0;
            value = va_arg(va, wchar_t*);
            int len = 0;
            for (; value[len] != '\0'; len++) {;}
            flags[4] = (flags[4] == -1) ? len : flags[4];
            for (int i = 0; i < flags[4] && value[i] != '\0'; i++, *index+=1) {
                buffer[*index] = value[i];
            }
        } else {
            char *value = 0;
            value = va_arg(va, char*);
            flags[4] = (flags[4] == -1) ? (int)s21_strlen(value) : flags[4];
            for (int i = 0; i < flags[4] && value[i] != '\0'; i++, *index+=1) {
                buffer[*index] = value[i];
            }
        }
    }
}
