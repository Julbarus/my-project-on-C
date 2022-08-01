#include "s21_string.h"

void processPointer(char* buffer, int* index, int flags[], void *value) {
    if (flags[0]) {
        *index = -1;
    } else if (flags[1] > 1 || flags[2]) {
        *index = -1;
    } else if (flags[4] != -1) {
        *index = -1;
    } else if ((flags[5] > 0 && flags[5] < 3) || flags[5] == 4) {
        *index = -1;
    } else if (*index != -1) {
        unsigned long adress = (unsigned long)value;
        char little_buff[12];
        int i = 0;
        buffer[*index] = '0';
        *index += 1;
        buffer[*index] = 'x';
        *index += 1;
        do {
            unsigned long uno = (int)(adress % 16);
            adress /= 16;
            little_buff[i] = uno > 9 ? (int)('W' + uno) : (int)('0' + uno);
            i++;
        } while (adress != 0);
        i--;
        for (; i >= 0; buffer[*index] = little_buff[i], *index += 1, i--) {;}
    }
}

void Space_16_8(char *buffer, const char *temp_buffer, int i, int *index, const int flags[]) {
    int width = flags[1] == 2 && flags[3] > i ? flags[3] - i - *index : flags[4] > i ? flags[4] - i : i;
    if (width > 0  && (flags[1] == 2 || flags[4] - i > 0)) {
        otstupSpase(buffer, index, width - 1, '0');
    }
    if (width > 0 || flags[4] != 0) {
        for (; i >= 0; buffer[*index] = temp_buffer[i], *index += 1, i--) {;}
    }
}

unsigned long int get_lengh(const int flags[], va_list va) {
    unsigned long int temp = flags[5] == 1 ? (unsigned short int)va_arg(va, int) :
    flags[5] == 3  || flags[5] ==4 ? va_arg(va, unsigned long int) :
    (unsigned long int)va_arg(va, unsigned int);
    return temp;
}

void processOctal(char* buffer, int* index, int flags[], va_list va) {
    char rBuffer[25] = "";
    int i = 0;
    unsigned long int n = get_lengh(flags, va);
    for (; n >= 8; i++, n /= 8) {
        int temp = n % 8;
        n -= temp;
        rBuffer[i] = temp + '0';
    }
    n < 8 ? rBuffer[i] = n + '0' : 0;
    flags[2] && n > 0 ? i++, rBuffer[i] = '0' : 0;
    Space_16_8(buffer, rBuffer, i, index, flags);
}

void processNum_16(char* buffer, int* index, int flags[], va_list va) {
    unsigned long int value = get_lengh(flags, va);
    char maxChr = flags[6] == 'x' ? 'W' : '7';
    char little_buff[254] = "";
    int i = 0;
    *index = flags[0] == 43 ? -1 : *index;
    for (; value >= 16 && *index != -1; i++) {
        unsigned uno = value % 16;
        value /= 16;
        little_buff[i] = uno > 9 ? maxChr + uno : '0' + uno;
    }
    if (value < 16 && *index != -1) {
        little_buff[i] = value > 9 ? maxChr + value : '0' + value;
    }
    if (flags[2] && value > 0 && *index != -1) {
        buffer[*index] = '0';
        *index+=1;
        buffer[*index] = flags[6] == 'x' ? 'x' : 'X';
        *index+=1;
    }
    *index != -1 ? Space_16_8(buffer, little_buff, i, index, flags) : 0;
}

void processN_Spec(int *value, const int *index) {
    *value = *index;
}
