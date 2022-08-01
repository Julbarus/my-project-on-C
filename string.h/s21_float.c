#include "s21_string.h"

char outputPlusOrMinus(long double *value, int* count) {
    char plusOrMinus;
    if (*value == 0) {
        *value = 0;
        plusOrMinus = '+';
    } else if (*value < 1.0) {
        do {
            *value *= 10;
            *count += 1;
        } while (*value < 0.9);
        plusOrMinus = '-';
    } else {
        while (*value >= 10.0) {
            *value /= (double)10;
            *count += 1;
        }
        plusOrMinus = '+';
    }
    return plusOrMinus;
}

long double processG(long double value, int flags[], int* count, char* plusOrMinus) {
    long double start_value = value;
    *plusOrMinus = outputPlusOrMinus(&value, count);
    if (flags[4] == 0)
        flags[4] = 1;
    if ((*count > 4 && *plusOrMinus == '-') ||
        (*count > 5 && *plusOrMinus == '+' && *count >= flags[4]) ||
        (amountNumber((int)start_value) > flags[4])) {
        if ((int)value != 0) {
            flags[4] -= amountNumber((int)value);
        }
    } else {
        long long int int_part = (long long int)start_value;
        long double fraction = start_value - (long double)int_part;
        if (int_part == 0) {
            for (int i = 0; i < 3; ++i) {
                if ((int)(fraction * 10) == 0) {
                    flags[4]++;
                    fraction *= 10;
                } else {
                    break;
                }
            }
        } else {
            flags[4] -= amountNumber(int_part);
        }
        flags[6] = 'F';
        value = start_value;
    }
    return value;
}

void outputFloat(char *buffer, long double fraction, int *index, int flags[]) {
    int int_part = fraction * 10;
    flags[4] += 1;
    fraction = (fraction * 10) - int_part;
    int inFlag = 0;
    for (int i = 0; (fraction > 0 && i < flags[4]) || i < flags[4];
        int_part = fraction * 10,
        fraction = (fraction * 10) - int_part, i++, *index += 1) {
        buffer[*index] = '0' + int_part;
        if (buffer[*index] > '4' && i == flags[4]-1) {
            int temp = *index-1;
            while (buffer[temp] == '9') {
                buffer[temp] = '0';
                temp--;
            }
            if (buffer[temp] =='.') {
                temp--;
                    while (buffer[temp] == '9') {
                    buffer[temp] = '0';
                    temp--;
                }
            }
            buffer[temp] = buffer[temp] + 1;
        }
        if (i == 1) inFlag = 1;
    }
    *index -= 1;
    if (inFlag && (flags[6] == 'g' || flags[6] == 'G' || flags[6] == 'F') && !flags[2]) {
        int i = *index - 1;
        for (; buffer[i] == '0' || buffer[i] == '.'; i--) {
            if (buffer[i] == '.') {
                buffer[i] = '\0';
                i--;
                break;
            }
            buffer[i] = '\0';
        }
        *index = i + 1;
    }
}

void expExp(char* tempBuff, int* tempIndex, const int flags[], char plusOrMinus, int count) {
    tempBuff[*tempIndex] = flags[6];
    *tempIndex+=1;
    tempBuff[*tempIndex] = plusOrMinus;
    *tempIndex+=1;
    if (count < 10) {
        tempBuff[*tempIndex] = '0';
        *tempIndex+=1;
    }
    for (int i = s21_powl(10, amountNumber(count) - 1); i > 0; i /= 10) {
        tempBuff[*tempIndex] = '0' + (count / i % 10);
        *tempIndex+=1;
    }
}

void zapolnenieFloatInt(char* tempBuff, int* tempIndex, int flags[], long double value) {
    int count = 0;
    char plusOrMinus;
    if (flags[4] == -1) flags[4] = 6;
    if (flags[6] == 'g' || flags[6] == 'G') {
        value = processG(value, flags, &count, &plusOrMinus);
    } else if (flags[6] == 'e' || flags[6] == 'E') {
        plusOrMinus = outputPlusOrMinus(&value, &count);
    }
    long long unsigned int int_part = (long long unsigned int)value;
    long double fraction = value - (long double)int_part;
    for (long long unsigned int i = s21_powl(10, amountNumber(int_part) - 1); i > 0; i /= 10) {
        tempBuff[*tempIndex] = '0' + (int_part / i % 10);
        *tempIndex+=1;
    }
    if (flags[4] > 0 || flags[2]) {
        tempBuff[*tempIndex] = '.';
        *tempIndex+=1;
    }
    outputFloat(tempBuff, fraction, tempIndex, flags);
    flags[6] == 'g' || flags[6] == 'G' ? flags[6] -= 2 : 0;
    if (flags[6] == 'e' || flags[6] == 'E') {
        expExp(tempBuff, tempIndex, flags, plusOrMinus, count);
    }
}

void nan_zapis(char *buffer, int *index) {
    buffer[*index] = 'n';
    *index += 1;
    buffer[*index] = 'a';
    *index += 1;
    buffer[*index] = 'n';
    *index += 1;
}

void inf_zapis(char *buffer, int *index) {
    buffer[*index] = 'i';
    *index += 1;
    buffer[*index] = 'n';
    *index += 1;
    buffer[*index] = 'f';
    *index += 1;
}

void processFloatAndExp(char* buffer, int* index, int flags[], va_list va) {
    long double value = 0;
    char tempBuff[1024];
    int tempIndex = 0;
    if (flags[5] == 2)
        value = va_arg(va, long double);
    else if (flags[5] == 1 || flags[5] == 4)
        *index = -1;
    else
        value = va_arg(va, double);
    int flag = isnan(value) ? 1 : isinf(value) ? 2 : 0;
    flag == 1 ? nan_zapis(buffer, index) : 0;
    flag == 2 ? inf_zapis(buffer, index) : 0;
    if (*index != -1 && flag == 0) {
        if (value < 0) {
            buffer[*index] = '-';
            *index+=1;
            value *= -1;
        } else if (flags[0]) {
            buffer[*index] = flags[0];
            *index+=1;
        }
        zapolnenieFloatInt(tempBuff, &tempIndex, flags, value);
            if (flags[1] == 2 && flags[3] > tempIndex) {
                otstupSpase(buffer, index, flags[3] - *index, '0');
                tempIndex < *index ? *index -= tempIndex : 0;
            }
        for (int i = 0; i < tempIndex; i++) {
            buffer[*index] = tempBuff[i];
            *index+=1;
        }
        buffer[*index] = '\0';
    }
}
