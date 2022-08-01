#include "s21_string.h"

void zapolnenie(char* buffer, int* index, int flags[], const char* tempBuff, const int *tempIndex) {
    if (flags[1] == 2 && flags[3] > *tempIndex) {
            otstupSpase(buffer, index, flags[3] - *index, '0');
        } else if (flags[4] > *tempIndex) {
            otstupSpase(buffer, index, flags[4], '0');
        }
        if (*tempIndex < *index) {
            *index-=*tempIndex;
        }
        for (int i = 0; i < *tempIndex; i++) {
            buffer[*index] = tempBuff[i];
            *index+=1;
        }
        buffer[*index] = '\0';
}

void zapolneniTempBufferaInt(char* buffer, int* index, int flags[], long long int value) {
    if (value == LLONG_MIN) {
        const char min[] = "-9223372036854775808";
        for (int i = 0; min[i] != '\0'; i++) {
            buffer[i] = min[i];
            *index+=1;
        }
    } else {
        if (value < 0) {
            buffer[*index] = '-';
            *index+=1;
            value *= -1;
        } else if (flags[0]) {
            buffer[*index] = flags[0];
            *index+=1;
        }
        char tempBuff[42] = "";
        int tempIndex = 0;
        for (long long int i = s21_powl(10, amountNumber(value) - 1); i > 0; i /= 10) {
            tempBuff[tempIndex] = '0' + (value / i % 10);
            tempIndex++;
        }
        zapolnenie(buffer, index, flags, tempBuff, &tempIndex);
    }
}

void zapolneniTempBufferaU(char* buffer, int* index, int flags[], long long unsigned int value) {
    char tempBuff[42] = "";
    int tempIndex = 0;
    if (flags[0]) {
        buffer[*index] = flags[0];
        *index+=1;
    }
    for (long long unsigned int i = s21_powl(10, amountNumber(value) - 1); i > 0; i /= 10) {
        tempBuff[tempIndex] = '0' + (value / i % 10);
        tempIndex++;
    }
    zapolnenie(buffer, index, flags, tempBuff, &tempIndex);
}

void processInt(char* buffer, int* index, int flags[], va_list va) {
    flags[2] ? *index = -1 : 0;
    if ((flags[6] == 'd' || flags[6] == 'i') && *index != -1) {
        long long int value = 0;
        if (flags[5] == 3) {
            value = va_arg(va, long int);
        } else if (flags[5] == 4 || (flags[5] == 2 && flags[6] != 'i')) {
            value = va_arg(va, long long int);
        } else if (flags[5] == 1) {
            value = (short)va_arg(va, int);
        } else {
            value = va_arg(va, int);
        }
        zapolneniTempBufferaInt(buffer, index, flags, value);
    }
}

void processUint(char* buffer, int* index, int flags[], va_list va) {
    flags[2] ? *index = -1 : 0;
    if (flags[6] == 'u' && *index != -1) {
        long long unsigned int value1 = 0;
        if (flags[5] == 1) {
            value1 = (unsigned short)va_arg(va, unsigned int);
        } else if (flags[5] == 2 || flags[5] == 4) {
            value1 = va_arg(va, long long unsigned int);
        } else if (flags[5] == 3) {
            value1 = va_arg(va, long unsigned int);
        } else {
            value1 = va_arg(va, unsigned int);
        }
        zapolneniTempBufferaU(buffer, index, flags, value1);
    }
}
