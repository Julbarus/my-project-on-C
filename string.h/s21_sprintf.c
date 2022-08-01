#include "s21_string.h"

long double s21_powl(long double a, int n) {
    long double res = 1.;
    while (n > 0) {
        if ((n & 1) != 0)
            res *= a;
        a *= a;
        n >>= 1;
    }
    return res;
}

int amountNumber(unsigned long long int n) {
    int count = 0;
    do {
        n /= 10;
        count++;
    } while (n != 0);
    return count;
}

void otstupSpase(char* buffer, int* index, int indexSpace, char c) {
    int sumIndex = indexSpace + *index;
    while (*index < sumIndex) {
        buffer[*index] = c;
        *index += 1;
    }
}

char* Accuracy_Place(char* ptr_format, int flags[], va_list va) {
    int ind = 0;
    int res = 0;
    if (*ptr_format == '0') {
        flags[1] = 2;
        ptr_format++;
    }
    if ((*ptr_format >= '0' && *ptr_format <= '9') || *ptr_format == '*') {
        ind = 3;
    } else if (*ptr_format == '.') {
        ind = 4;
        ptr_format++;
        flags[ind] = 0;
    }
    if (*ptr_format == '*') {
        ptr_format++;
        flags[ind] = va_arg(va, int);
        if (flags[ind] < 0) {
            flags[1] = 1;
            flags[ind] *= -1;
        }
    } else if (*ptr_format >= '0' && *ptr_format <= '9') {
        while (*ptr_format >= '0' && *ptr_format <= '9') {
            res = res * 10 + (*ptr_format - 48);
            ptr_format++;
        }
        flags[ind] = res;
    }
    return ptr_format;
}

char* findFlags(int flags[], int* indexForBuffer, char* ptr_format, va_list va) {  // До этой нужны
    while (*ptr_format && *indexForBuffer != -1) {
        if ((*ptr_format == ' ' || *ptr_format == '+') && !flags[0]) {
            flags[0] = *ptr_format;
            ptr_format++;
        } else if (*ptr_format == '-' && !flags[1]) {
            flags[1] = 1;
            ptr_format++;
        } else if (*ptr_format == '#') {
            flags[2] = 1;
            ptr_format++;
        } else if (*ptr_format == 'h' && !flags[5]) {
            flags[5] = 1;
            ptr_format++;
        } else if (*ptr_format == 'L' && !flags[5]) {
            flags[5] = 2;
            ptr_format++;
        } else if (*ptr_format == 'l') {
            if (flags[5] == 3) {
                flags[5] = 4;
            } else if (flags[5] == 0) {
                flags[5] = 3;
            } else {
                *indexForBuffer = -1;
            }
            ptr_format++;
        } else if ((*ptr_format >= '0' && *ptr_format <= '9') || (*ptr_format == '*' || *ptr_format == '.')) {
            ptr_format = Accuracy_Place(ptr_format, flags, va);
        } else if (s21_strchr("difcsueEgGxXonp%", *ptr_format) != s21_NULL) {
            flags[6] = *ptr_format++;
            break;
        } else {
            *indexForBuffer = -1;
        }
    }
    return ptr_format;
}


void specificator(char* buffer, int* indexForBuffer, int flags[], va_list va) {  // Обработка конфликтов тут!
    if (flags[6] == 'd' || flags[6] == 'i') {
        processInt(buffer, indexForBuffer, flags, va);
    } else if (flags[6] == 'u') {
        processUint(buffer, indexForBuffer, flags, va);
    } else if (flags[6] == 'f' || flags[6] == 'e' || flags[6] == 'E' || flags[6] == 'g' || flags[6] == 'G') {
        processFloatAndExp(buffer, indexForBuffer, flags, va);
    } else if (flags[6] == 'c') {
        processChar(buffer, indexForBuffer, flags, va);
    } else if (flags[6] == 's') {
        processStr(buffer, indexForBuffer, flags, va);
    } else if (flags[6] == 'o') {
        processOctal(buffer, indexForBuffer, flags, va);
    } else if (flags[6] == 'p') {
        processPointer(buffer, indexForBuffer, flags, va_arg(va, void*));
    } else if (flags[6] == '%') {
        for (int i = 0; i < 6; i++) {
            flags[i] > 0 ? *indexForBuffer = -1 : 0;
        }
        buffer[*indexForBuffer] = '%';
        *indexForBuffer+=1;
    } else if (flags[6] == 'x' || flags[6] == 'X') {
        processNum_16(buffer, indexForBuffer, flags, va);
    } else {
        *indexForBuffer = -1;
    }
}

int parser(char* buffer, char* format, va_list va) {
    int index_for_buffer = 0;
    char *ptr_format = format;
    // i = 0 ' ' = 32 and '+' = 43
    // i = 1 '-' = 1 and '0' = 2 если '0' то он смотрит на 3 ячейку
    // i = 2 '#'
    // i = 3 'число and *'
    // i = 4 '.число and .*'
    // i = 5 'h = 1 AND L = 2 and l = 3 and ll = 4'
    // i = 6 'Спецификатор'
    // '0', '#' и '+' Только для числовых переменных
    // Точность только для числовых и строки
    // h, L и l тоже только для числовых переменных
    while (*ptr_format && index_for_buffer != -1) {
        char tempBuffer[1024] = {0};
        int tempIndexBuffer = 0;
        if (*ptr_format != '%') {
            buffer[index_for_buffer++] = *ptr_format++;
            continue;
        }
        int flags[] = {0, 0, 0, 0, -1, 0, 0};
        ptr_format++;
        ptr_format = findFlags(flags, &tempIndexBuffer, ptr_format, va);
        if (flags[6] == 'n') {
            processN_Spec(va_arg(va, int *), &index_for_buffer);
            continue;
        }
        specificator(tempBuffer, &tempIndexBuffer, flags, va);
        if (!flags[1] && flags[3] > tempIndexBuffer && tempIndexBuffer != -1) {
            otstupSpase(buffer, &index_for_buffer, flags[3], ' ');
            index_for_buffer -= tempIndexBuffer;
        }
        for (int i = 0; i < tempIndexBuffer; i++) {
            buffer[index_for_buffer] = tempBuffer[i];
            index_for_buffer++;
        }
        if (flags[1] == 1 && flags[3] > tempIndexBuffer) {
            otstupSpase(buffer, &index_for_buffer, flags[3] - tempIndexBuffer, ' ');
        }
        if (tempIndexBuffer == -1) {
            index_for_buffer = tempIndexBuffer;
            buffer[0] = '\0';
        }
    }
    if (index_for_buffer > 0) {
        buffer[index_for_buffer] = '\0';
    }
    return index_for_buffer;
}

int s21_sprintf(char* buffer, char* format, ...) {
    va_list va;
    va_start(va, format);
    const int ret = parser(buffer, format, va);
    va_end(va);
    return ret;
}
