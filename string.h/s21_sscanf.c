#include "s21_string.h"

const char* processCharScn(int* amount, const char* ptr_str, const int flags[], va_list va, int* flagCount) {
    if (flags[2] && *ptr_str != '\0') {
        if (flags[0] == 3)
            *va_arg(va, wchar_t *) = *ptr_str;
        else
            *va_arg(va, char *) = *ptr_str;
        *amount = *flagCount ? *amount + 1 : *amount;
    }
    ptr_str++;
    return ptr_str;
}

int strValidScn(char c) {
    int res = 1;
    switch (c) {
        case '\0':
        case '\n':
        case '\t':
        case '\v':
        case '\f':
        case '\r':
        case ' ':
            res = 0;
            break;
    }
    return res;
}

const char* getStrWchar(const char* ptr_str, int flags[], va_list va) {
    wchar_t *temp_ptr = va_arg(va, wchar_t *);
    while (strValidScn(*ptr_str) && flags[1] > 0) {
        *temp_ptr = *ptr_str;
        temp_ptr++;
        ptr_str++;
        flags[1]--;
    }
    *temp_ptr = '\0';
    return ptr_str;
}

const char* getStr(const char* ptr_str, int flags[], va_list va) {
    char *temp_ptr = va_arg(va, char *);
    while (strValidScn(*ptr_str) && flags[1] > 0) {
        *temp_ptr = *ptr_str;
        temp_ptr++;
        ptr_str++;
        flags[1]--;
    }
    *temp_ptr = '\0';
    return ptr_str;
}

const char* processStrScn(int* amount, const char* ptr_str, int flags[], va_list va, int* flagCount) {
    if (flags[2] && strValidScn(*ptr_str)) {
        if (flags[0] == 3)
            ptr_str = getStrWchar(ptr_str, flags, va);
        else
            ptr_str = getStr(ptr_str, flags, va);
        *amount = *flagCount ? *amount + 1 : *amount;
    } else {
        while (*ptr_str != ' ')
            ptr_str++;
    }
    return ptr_str;
}

int convertScn(int value) {
    if (value >= '0' && value <= '9')
        value -= '0';
    else if (value == 'x')
        value = 0;
    else if (value == 'a' || value == 'A')
        value = 10;
    else if (value == 'b' || value == 'B')
        value = 11;
    else if (value == 'c' || value == 'C')
        value = 12;
    else if (value == 'd' || value == 'D')
        value = 13;
    else if (value == 'e' || value == 'E')
        value = 14;
    else if (value == 'f' || value == 'F')
        value = 15;
    return value;
}

int validScn(int num, const char* ptr_str) {
    int res = 0;
    if (num == 10 && (*ptr_str >= '0' && *ptr_str <= '9')) res = 1;
    if (num == 8 && (*ptr_str >= '0' && *ptr_str <= '7')) res = 1;
    if (num == 16 && ((*ptr_str >= '0' && *ptr_str <= '9') ||
        (*ptr_str >= 'a' && *ptr_str <= 'f') ||
        (*ptr_str >= 'A' && *ptr_str <= 'F') ||
        (*ptr_str == 'x' || *ptr_str == 'X'))) res = 1;
    return res;
}

const char* ifMinus(int flags[], const char* ptr_str) {
    if (*ptr_str == '-') {
        flags[4] = 1;
        flags[1]--;
        ptr_str++;
    }
    return ptr_str;
}

void countVariable(int count, int *flagCount, const int flags[], int *amount) {
    if (count && *flagCount && flags[2])
        *amount += 1;
    else if (!count)
        *flagCount = 0;
}

void runProcess16And8(int flags[], const char* ptr_str, int *count, int *flagIn) {
    if ((*count == 0 && *ptr_str == '0') || (*count == 0 && (*ptr_str == 'x' || *ptr_str == 'X')))
        *count = 0;
    else
        *count += 1;
    *flagIn = 1;
    flags[1]--;
}

const char* getValueNumScn(int flags[], const char* ptr_str, int* amount, int num,
    int* flagCount, va_list va) {
    unsigned int value = 0;
    ptr_str = ifMinus(flags, ptr_str);
    int flagIn = 0;
    int count = 0;
    if (validScn(num, ptr_str)) {
        while (validScn(num, ptr_str) && flags[1] > 0) {
            runProcess16And8(flags, ptr_str, &count, &flagIn);
            value = value * num + convertScn(*ptr_str);
            ptr_str++;
        }
    } else {
        flags[2] = 0;
    }
    if ((flags[3] == 'x' || flags[3] == 'X') && count > 16) value = -1;
    countVariable(flagIn, flagCount, flags, amount);
    if (flags[4]) value *= -1;
    if (flags[2] && flags[0] == 0)
        *va_arg(va, unsigned int *) = value;
    else if (flags[2] && flags[0] == 1)
        *va_arg(va, unsigned short int *) = (unsigned short int)value;
    return ptr_str;
}

const char* getValueLongNumScn(int flags[], const char* ptr_str, int* amount, int num,
    int* flagCount, va_list va) {
    unsigned short value = 0;
    ptr_str = ifMinus(flags, ptr_str);
    int flagIn = 0;
    int count = 0;
    if (validScn(num, ptr_str)) {
        while (validScn(num, ptr_str) && flags[1] > 0) {
            runProcess16And8(flags, ptr_str, &count, &flagIn);
            value = value * num + convertScn(*ptr_str);
            ptr_str++;
        }
    } else {
        flags[2] = 0;
    }
    if ((flags[3] == 'x' || flags[3] == 'X') && count > 16) value = -1;
    countVariable(flagIn, flagCount, flags, amount);
    if (flags[4]) value *= -1;
    if (flags[2]) *va_arg(va, unsigned long int *) = value;
    return ptr_str;
}

const char* processNum_8_or_16Scn(int* amount, const char* ptr_str, int flags[], va_list va, int* flagCount) {
    if (flags[3] == 'o') {
        if (flags[0] == 1 || flags[0] == 0)
            ptr_str = getValueNumScn(flags, ptr_str, amount, 8, flagCount, va);
        else
            ptr_str = getValueLongNumScn(flags, ptr_str, amount, 8, flagCount, va);
    } else {
        if (flags[0] == 1 || flags[0] == 0)
            ptr_str = getValueNumScn(flags, ptr_str, amount, 16, flagCount, va);
        else
            ptr_str = getValueLongNumScn(flags, ptr_str, amount, 16, flagCount, va);
    }
    return ptr_str;
}

const char* getValueLongLongScn(const char* ptr_str, int flags[], long long* value,
    int* amount, int* flagCount) {
    long long unsigned tempUn = 0;
    ptr_str = ifMinus(flags, ptr_str);
    int count = 0;
    if (validScn(10, ptr_str)) {
        while (validScn(10, ptr_str) && flags[1] > 0) {
            count = 1;
            tempUn = tempUn * 10 + (*ptr_str - 48);
            ptr_str++;
            flags[1]--;
            if (flags[3] == 'd' && tempUn > LLONG_MAX) {
                tempUn = LLONG_MAX;
                break;
            }
        }
    } else {
        flags[2] = 0;
    }
    *value = tempUn;
    countVariable(count, flagCount, flags, amount);
    if (flags[4]) *value *= -1;
    return ptr_str;
}

const char* processD_and_UScn(int* amount, const char* ptr_str, int flags[], va_list va, int* flagCount) {
    long long value = 0;
    ptr_str = getValueLongLongScn(ptr_str, flags, &value, amount, flagCount);
    if (flags[2]) {
        if (flags[3] == 'd') {
            if (flags[0] == 3) {
                *va_arg(va, long int *) = (long)value;
            } else if (flags[0] == 4 || (flags[0] == 2)) {
                *va_arg(va, long long int *) = value;
            } else if (flags[0] == 1) {
                *va_arg(va, short *) = (short)value;
            } else {
                *va_arg(va, int *) = (int)value;
            }
        } else if (flags[3] == 'u') {
            if (flags[0] == 3) {
                *va_arg(va, long unsigned *) = (long unsigned)value;
            } else if (flags[0] == 4 || (flags[0] == 2)) {
                *va_arg(va, long long unsigned *) = value;
            } else if (flags[0] == 1) {
                *va_arg(va, unsigned short *) = (unsigned short)value;
            } else {
                *va_arg(va, unsigned int *) = (unsigned int)value;
            }
        }
    }
    return ptr_str;
}

const char* processIScn(int* amount, const char* ptr_str, int flags[], va_list va, int* flagCount) {
    ptr_str = ifMinus(flags, ptr_str);
    if (*ptr_str == '0') {
        ptr_str++;
        if (*ptr_str == 'x' || *ptr_str == 'X') {
            flags[3] = 'x';
        } else {
            ptr_str--;
            flags[3] = 'o';
        }
        ptr_str = processNum_8_or_16Scn(amount, ptr_str, flags, va, flagCount);
    } else {
        flags[3] = 'd';
        ptr_str = processD_and_UScn(amount, ptr_str, flags, va, flagCount);
    }
    return ptr_str;
}

const char* getFractionPartScn(const char* ptr_str, int flags[], long double* value) {
    int count = 0;
    while ((*ptr_str >= '0' && *ptr_str <= '9') && flags[1] > 0) {
        count++;
        *value = *value + (*ptr_str - 48) / pow(10, count);
        ptr_str++;
        flags[1]--;
    }
    return ptr_str;
}

const char* getFloat(int* amount, const char* ptr_str, int flags[], int* flagCount,
    long long int *int_part, long double *fraction, int *exp, int *plusOrMinus) {
    ptr_str = getValueLongLongScn(ptr_str, flags, int_part, amount, flagCount);
    if (flags[4]) *int_part *= -1;
    if (*ptr_str == '.') {
        flags[1]--;
        ptr_str++;
        ptr_str = getFractionPartScn(ptr_str, flags, fraction);
    }
    if ((*ptr_str == 'e' || *ptr_str == 'E') && flags[1] > 1) {
        int flagNumber = 0, flagPlus = 0;
        flags[1]--;
        ptr_str++;
        flagNumber++;
        if (((*ptr_str == '-' || *ptr_str == '+') || (*ptr_str >= '0' && *ptr_str <= '9'))
            && flags[1] > 0) {
            if (*ptr_str == '-' || *ptr_str == '+') {
                *plusOrMinus = (*ptr_str == '-') ? 0 : 1;
                ptr_str++;
                flagPlus++;
                flagNumber++;
                flags[1]--;
            }
            *plusOrMinus = (*plusOrMinus == 0) ? 0 : 1;
            while ((*ptr_str >= '0' && *ptr_str <= '9') && flags[1] > 0) {
                if (flagPlus) flagNumber++;
                *exp = *exp * 10 + (*ptr_str - 48);
                ptr_str++;
                flags[1]--;
            }
        }
        if (flagNumber == 2 && flagPlus) {
            for (int i = 0; i < 2; i++)
                ptr_str--;
        }
    }
    return ptr_str;
}

const char* processFloatAndExpScn(int* amount, const char* ptr_str, int flags[], va_list va, int* flagCount) {
    long long int int_part = 0;
    long double fraction = 0;
    int exp = 0, plusOrMinus = -1;  // 1 + 0 -
    ptr_str = getFloat(amount, ptr_str, flags, flagCount, &int_part, &fraction, &exp, &plusOrMinus);
    double value = (int_part + fraction);
    if (plusOrMinus != -1) {
        exp = (plusOrMinus == 1) ? exp : exp * (-1);
        value *= pow(10, exp);
    }
    if (flags[4]) value *= -1;
    if (flags[2]) {
        if (flags[0] == 3) {
            *va_arg(va, double *) = (double)value;
        } else if (flags[0] == 2) {
            *va_arg(va, long double *) = (long double)value;
        } else {
            *va_arg(va, float *) = (float)value;
        }
    }
    return ptr_str;
}

const char* processPointerScn(int* amount, const char* ptr_str, int flags[], va_list va, int* flagCount) {
    void **path = va_arg(va, void**);
    long adres = 0;
    int point = 1;
    if (validScn(16, ptr_str)) {
        while (validScn(16, ptr_str) && flags[1] > 0) {
            point = 0;
            adres = adres * 16 + convertScn(*ptr_str);
            ptr_str++;
            flags[1]--;
        }
    } else {
        flags[2] = 0;
    }
    if (flags[2]) {
        if (point)
            *path = &path;
        else
            *path = (void*)adres;
        *amount = *flagCount ? *amount + 1 : *amount;
    }
    return ptr_str;
}

const char* processProc(const char* ptr_str, int* flagCount) {
    if (*ptr_str == '%') {
        ptr_str++;
    } else {
        *flagCount = 0;
    }
    return ptr_str;
}

const char* specificatorScn(int* amount, const char* ptr_str, int flags[], va_list va, int *flagCount) {
    if (flags[3] == 'd') {
        ptr_str = processD_and_UScn(amount, ptr_str, flags, va, flagCount);
    } else if (flags[3] == 'i') {
        ptr_str = processIScn(amount, ptr_str, flags, va, flagCount);
    } else if (flags[3] == 'u') {
        ptr_str = processD_and_UScn(amount, ptr_str, flags, va, flagCount);
    } else if (flags[3] == 'f' || flags[3] == 'e' || flags[3] == 'E' ||
        flags[3] == 'g' || flags[3] == 'G') {
        ptr_str =  processFloatAndExpScn(amount, ptr_str, flags, va, flagCount);
    } else if (flags[3] == 'c') {
        ptr_str = processCharScn(amount, ptr_str, flags, va, flagCount);
    } else if (flags[3] == 's') {
        ptr_str = processStrScn(amount, ptr_str, flags, va, flagCount);
    } else if (flags[3] == 'o') {
        ptr_str = processNum_8_or_16Scn(amount, ptr_str, flags, va, flagCount);
    } else if (flags[3] == 'p') {
        ptr_str = processPointerScn(amount, ptr_str, flags, va, flagCount);
    } else if (flags[3] == '%') {
        ptr_str = processProc(ptr_str, flagCount);
    } else if (flags[3] == 'x' || flags[3] == 'X') {
        ptr_str = processNum_8_or_16Scn(amount, ptr_str, flags, va, flagCount);
    }
    return ptr_str;
}

const char* widthScn(const char* ptr_format, int flags[]) {
    int res = 0;
    if (*ptr_format == '*') {
        ptr_format++;
        flags[2] = 0;
    } else if (*ptr_format >= '0' && *ptr_format <= '9') {
        while (*ptr_format >= '0' && *ptr_format <= '9') {
            res = res * 10 + (*ptr_format - 48);
            ptr_format++;
        }
        flags[1] = res;
    }
    return ptr_format;
}

const char* findFlagsScn(int flags[], const char* ptr_format) {
    while (*ptr_format) {
        if (*ptr_format == 'h' && !flags[0]) {
            flags[0] = 1;
            ptr_format++;
        } else if (*ptr_format == 'L' && !flags[0]) {
            flags[0] = 2;
            ptr_format++;
        } else if (*ptr_format == 'l') {
            if (flags[0] == 3) {
                flags[0] = 4;
            } else if (flags[0] == 0) {
                flags[0] = 3;
            }
            ptr_format++;
        } else if ((*ptr_format >= '0' && *ptr_format <= '9') || (*ptr_format == '*')) {
            ptr_format = widthScn(ptr_format, flags);
        } else if (s21_strchr("difcsueEgGxXonp%", *ptr_format) != NULL) {
            flags[3] = *ptr_format++;
            break;
        }
    }
    if (!flags[1]) flags[1] = 512;
    return ptr_format;
}

// i = 0 lenght 'h = 1 AND L = 2 and l = 3 and ll = 4'
// i = 1 widthScn
// i = 2 write flag
// i = 3 specificatorScn
// i = 4 minus
int parserScn(const char* str, const char* format, va_list va) {
    const char *ptr_str = str;
    const char *ptr_format = format;
    int amount = *ptr_str ? 0 : -1;
    int flagCount = 1;
    while (*ptr_format) {
        if (*ptr_format == '%') {
            int flags[] = {0, 0, 1, 0, 0};
            ptr_format++;
            ptr_format = findFlagsScn(flags, ptr_format);
            if (flags[3] == 'n') {
                *va_arg(va, int *) = (int)(ptr_str - str);
            } else {
                if (flags[3] != 'c')
                    while (*ptr_str == ' ')
                        ptr_str++;
                ptr_str = specificatorScn(&amount, ptr_str, flags, va, &flagCount);
            }
        } else if (*ptr_format == ' ') {
            while (*ptr_format == ' ')
                ptr_format++;
            while (*ptr_str == ' ')
                ptr_str++;
        }
        if (!flagCount)
            break;
    }
    return amount;
}

int s21_sscanf(const char* str, const char* format, ...) {
    va_list va;
    va_start(va, format);
    int ret = parserScn(str, format, va);
    va_end(va);
    return ret;
}
