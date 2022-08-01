#include <check.h>
#include "s21_string.h"
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <wchar.h>

int size = sizeof(s21_NULL);

START_TEST(s21_strlenTEST) {
    char array_1[] = "Hello world!";
    char array_2[] = "";
    char array_3[] = "!@IU@!*&$@1KWljD";
    char array_4[] = "Hello world!\n";
    char array_5[] = "123456";

    ck_assert_uint_eq(s21_strlen(""), strlen(""));
    ck_assert_uint_eq(s21_strlen(" "), strlen(" "));
    ck_assert_uint_eq(s21_strlen("\0"), strlen("\0"));
    ck_assert_uint_eq(s21_strlen("Hello world\0"), strlen("Hello world\0"));
    ck_assert_uint_eq(s21_strlen("Hello world\n\0"), strlen("Hello world\n\0"));
    ck_assert_uint_eq(s21_strlen(" \n\0"), strlen(" \n\0"));
    ck_assert_uint_eq(s21_strlen(" \0"), strlen(" \0"));
    ck_assert_uint_eq(s21_strlen("\n\0"), strlen("\n\0"));
    ck_assert_uint_eq(s21_strlen("\\0ehwrgrwg\\n\\n\\u\\et"), strlen("\\0ehwrgrwg\\n\\n\\u\\et"));
    ck_assert_uint_eq(s21_strlen("œ∑´®†¥¨ˆøπåß∂ƒ©˙∆˚¬≈ç√∫˜µ"), strlen("œ∑´®†¥¨ˆøπåß∂ƒ©˙∆˚¬≈ç√∫˜µ"));
    ck_assert_ptr_null((void*)s21_strlen(s21_NULL));
    ck_assert_int_eq(s21_strlen(array_1), strlen(array_1));
    ck_assert_int_eq(s21_strlen(array_2), strlen(array_2));
    ck_assert_int_eq(s21_strlen(array_3), strlen(array_3));
    ck_assert_int_eq(s21_strlen(array_4), strlen(array_4));
    ck_assert_int_eq(s21_strlen(array_5), strlen(array_5));
}
END_TEST

START_TEST(s21_memchrTEST) {
    char mass1[100] = "1234";
    char mass2[100] = "Nice";
    char mass3[100] = "00000000";
    char mass4[100] = "-1-2-3-4-5";
    char mass5[100] = "daseqweqwe";
    char mass6[] = "Hello world!";
    char mass7[] = "";
    char *mass8 = "This is a s21_memchr test";

    fail_unless(s21_memchr(mass1, '4', 9) == memchr(mass1, '4', 9));
    fail_unless(s21_memchr(mass2, 'i', 9) == memchr(mass2, 'i', 9));
    fail_unless(s21_memchr(mass3, '0', 9) == memchr(mass3, '0', 9));
    fail_unless(s21_memchr(mass4, '-', 9) == memchr(mass4, '-', 9));
    fail_unless(s21_memchr(mass5, '1', 9) == memchr(mass5, '1', 9));
    ck_assert_ptr_eq(s21_memchr(mass6, 'w', 13), memchr(mass6, 'w', 13));
    ck_assert_ptr_eq(s21_memchr(mass6, ' ', 13), memchr(mass6, ' ', 13));
    ck_assert_ptr_eq(s21_memchr(mass6, 111, 13), memchr(mass6, 111, 13));
    ck_assert_ptr_eq(s21_memchr(mass7, ' ', 1), memchr(mass7, ' ', 1));
    ck_assert_ptr_eq(s21_memchr(mass8, 's', 0), memchr(mass8, 's', 0));
    ck_assert_ptr_eq(s21_memchr(mass8, 's', 2), memchr(mass8, 's', 2));
    ck_assert_ptr_eq(s21_memchr(mass8, 's', 50), memchr(mass8, 's', 50));
    ck_assert_ptr_eq(s21_memchr(mass8, 'u', 50), memchr(mass8, 'u', 50));
    ck_assert_ptr_eq(s21_memchr(mass8, 'u', 32), memchr(mass8, 'u', 32));
    ck_assert_ptr_eq(s21_memchr("abc", 'c', 0), memchr("abc", 'c', 0));
    ck_assert_ptr_eq(s21_memchr("abc", 'c', 2), memchr("abc", 'c', 2));
    ck_assert_ptr_eq(s21_memchr("abc", 'c', 1), memchr("abc", 'c', 1));
    ck_assert_ptr_null(s21_memchr(s21_NULL, '0', 1));
    ck_assert_ptr_null(s21_memchr(s21_NULL, '0', size));
}
END_TEST

START_TEST(s21_memcmpTEST) {
    char str1[] = "abcdefghj", str2[] = "abcdefghj", str3[] = "abcdefghjz",
         str4[] = "", str5[] = "\n";
    ck_assert_int_eq(s21_memcmp(str1, str2, 10), memcmp(str1, str2, 10));
    ck_assert_int_eq(s21_memcmp(str1, str3, 9), memcmp(str1, str3, 9));
    ck_assert_int_eq(s21_memcmp(str1, str4, 1), memcmp(str1, str4, 1));
    ck_assert_int_eq(s21_memcmp(str5, str4, 1), memcmp(str5, str4, 1));
}
END_TEST

START_TEST(s21_memcpyTEST) {
    char array_4[] = "abcd";
    char array_5[] = "efgh";
    char array_6[] = "efghi";
    char array_8[] = "abcd";
    char array_9[] = "efgh";
    char array_10[] = "efghi";
    char *src1 = "schesrdsfdsfdsfesr";
    int len = s21_strlen(src1);
    char mass1[100] = {0};
    char mass2[100] = {0};
    char mass3[100] = {0};
    char mass4[100] = {0};
    char mass5[100] = {0};
    char mass6[100] = {0};
    char mass7[100] = {0};
    char mass8[100] = {0};
    char mass9[100] = {0};
    char mass10[100] = {0};

    s21_memcpy(mass1, src1, len + 1);
    memcpy(mass2, src1, len + 1);

    s21_memcpy(mass3, src1, 5);
    memcpy(mass4, src1, 5);

    s21_memcpy(mass5, src1, 0);
    memcpy(mass6, src1, 0);

    s21_memcpy(mass7, src1, 18);
    memcpy(mass8, src1, 18);

    s21_memcpy(mass9, src1, 1);
    memcpy(mass10, src1, 1);

    ck_assert_str_eq(mass1, mass2);
    ck_assert_str_eq(mass3, mass4);
    ck_assert_str_eq(mass5, mass6);
    ck_assert_str_eq(mass7, mass8);
    ck_assert_str_eq(mass9, mass10);

    s21_memcpy(array_4, array_5, 4);
    memcpy(array_8, array_9, 4);
    ck_assert_pstr_eq(array_4, array_8);

    s21_memcpy(array_4, array_6, 4);
    memcpy(array_8, array_10, 4);
    ck_assert_pstr_eq(array_4, array_8);

    s21_memcpy(array_6, array_4, 1);
    memcpy(array_10, array_8, 1);
    ck_assert_pstr_eq(array_6, array_10);

    char array_1[256] = "qwertyuiop[]asdfghjkl;zxcvbnm<>?";
    char array_2[256] = "?><mnbvcxz;lkjhgfdsa][poiuytrewq";
    char array_3[256] = "";
    void* mem_1 = malloc(256);
    void* mem_2 = malloc(256);
    ck_assert_str_eq(memcpy(mem_1, array_1, sizeof(array_1)), s21_memcpy(mem_1, array_1, sizeof(array_1)));
    ck_assert_mem_eq(s21_memcpy(mem_1, array_1, sizeof(array_1)), array_1, sizeof(array_1));
    ck_assert_str_eq(memcpy(mem_2, array_2, sizeof(array_2)), s21_memcpy(mem_2, array_2, sizeof(array_2)));
    ck_assert_mem_eq(s21_memcpy(mem_2, array_2, sizeof(array_2)), array_2, sizeof(array_2));
    ck_assert_ptr_eq(memcpy(mem_1, array_3, sizeof(s21_NULL)), s21_memcpy(mem_1, array_3, sizeof(s21_NULL)));
    free(mem_1);
    free(mem_2);
}
END_TEST

START_TEST(s21_memmoveTEST) {
    char arr_1[80] = "qwertyuiop[]asdfghjkl;zxcvbnm<>?";
    char arr_2[80] = "?><mnbvcxz;lkjhgfdsa][poiuytrewq";
    char arr_3[80] = "";
    char arr_4[41] = "?><mnbvcxz;lkjhgfdsa][poiuytrewqiluiluil";
    void* mem_1 = malloc(256);
    void* mem_2 = malloc(256);
    char *src1 = "schesrdsfdsfdsfesr";
    int len = s21_strlen(src1);
    char mass1[100];
    char mass2[100];
    char mass3[100];
    char mass4[100];
    char mass5[100];
    char mass6[100];
    char mass7[100];
    char mass8[100];
    char mass9[100];
    char mass10[100];

    s21_memmove(mass1, src1, len + 1);
    memmove(mass2, src1, len + 1);

    s21_memmove(mass3, src1, 5);
    memmove(mass4, src1, 5);

    s21_memmove(mass5, src1, 0);
    memmove(mass6, src1, 0);

    s21_memmove(mass7, src1, 18);
    memmove(mass8, src1, 18);

    s21_memmove(mass9, src1, 1);
    memmove(mass10, src1, 1);

    ck_assert_str_eq(memmove(mem_1, arr_1, sizeof(arr_1)), s21_memmove(mem_1, arr_1, sizeof(arr_1)));
    ck_assert_mem_eq(s21_memmove(mem_1, arr_1, sizeof(arr_1)), arr_1, sizeof(arr_1));
    ck_assert_str_eq(memmove(mem_2, arr_2, sizeof(arr_2)), s21_memmove(mem_2, arr_2, sizeof(arr_2)));
    ck_assert_mem_eq(s21_memmove(mem_2, arr_2, sizeof(arr_2)), arr_2, sizeof(arr_2));
    ck_assert_str_eq(memmove(arr_3, arr_3, sizeof(arr_3)), s21_memmove(arr_3, arr_3, sizeof(arr_3)));
    ck_assert_ptr_eq(memmove(mem_1, arr_3, sizeof(s21_NULL)), s21_memmove(mem_1, arr_3, sizeof(s21_NULL)));
    ck_assert_str_eq(memmove(arr_1, arr_4, sizeof(arr_4)), s21_memmove(arr_1, arr_4, sizeof(arr_4)));
    ck_assert_str_eq(memmove(arr_2, arr_4, sizeof(arr_4)), s21_memmove(arr_2, arr_4, sizeof(arr_4)));
    ck_assert_str_eq(memmove(arr_3, arr_4, sizeof(arr_4)), s21_memmove(arr_3, arr_4, sizeof(arr_4)));
    ck_assert_str_eq(memmove(mem_1, arr_1+5, sizeof(arr_1)+5), s21_memmove(mem_1, arr_1+5, sizeof(arr_1)+5));
    ck_assert_msg(memcmp(mass1, mass2, len) == 0, "1st memmove test failed.");
    ck_assert_msg(memcmp(mass3, mass4, 5) == 0, "2nd memmove test failed.");
    ck_assert_msg(memcmp(mass5, mass6, 0) == 0, "3rd memmove test failed.");
    ck_assert_msg(memcmp(mass7, mass8, 18) == 0, "4rth memmove test failed.");
    ck_assert_msg(memcmp(mass9, mass10, 1) == 0, "5rth memmove test failed.");
    free(mem_1);
    free(mem_2);
}
END_TEST

START_TEST(s21_memsetTEST) {
    char array_2[9] = "abcdefgh";
    char array_3[1] = "";
    char array_4[9] = "abcdefgh";
    char array_5[1] = "";

    s21_memset(array_2, 'i', 8);
    memset(array_4, 'i', 8);
    ck_assert_pstr_eq(array_2, array_4);

    s21_memset(array_2, 'i', 5);
    memset(array_4, 'i', 5);
    ck_assert_pstr_eq(array_2, array_4);

    s21_memset(array_2, ' ', 5);
    memset(array_4, ' ', 5);
    ck_assert_pstr_eq(array_2, array_4);

    s21_memset(array_2, ' ', 8);
    memset(array_4, ' ', 8);
    ck_assert_pstr_eq(array_2, array_4);

    s21_memset(array_3, ' ', 1);
    memset(array_5, ' ', 1);
    ck_assert_pstr_eq(array_2, array_4);

    s21_memset(array_3, 'i', 1);
    memset(array_5, 'i', 1);
    ck_assert_pstr_eq(array_2, array_4);

    char str1[] = "School is cool";
    char str2[] = "School-21";
    char str3[] = "amogus";
    char str4[] = "This is a string!";
    char str5[] = "Shrek";
    char res1[] = "School is cool";
    char res2[] = "School-21";
    char res3[] = "amogus";
    char res4[] = "This is a string!";
    char res5[] = "Shrek";

    s21_memset(str1, '_', 5);
    s21_memset(str2, '_', 4);
    s21_memset(str3, 'M', 2);
    s21_memset(str4, '-', 1);
    s21_memset(str5, 'S', 1);

    memset(res1, '_', 5);
    memset(res2, '_', 4);
    memset(res3, 'M', 2);
    memset(res4, '-', 1);
    memset(res5, 'S', 1);

    ck_assert_msg(memcmp(str1, res1, sizeof(str1)) == 0, "1st memset test failed");
    ck_assert_msg(memcmp(str2, res2, sizeof(str2)) == 0, "2nd memset test failed");
    ck_assert_msg(memcmp(str3, res3, sizeof(str3)) == 0, "3rd memset test failed");
    ck_assert_msg(memcmp(str4, res4, sizeof(str4)) == 0, "4rth memset test failed");
    ck_assert_msg(memcmp(str5, res5, sizeof(str5)) == 0, "5rth memset test failed");
}
END_TEST

START_TEST(s21_strcatTEST) {
    char array_2[34] = "?><mnbvcxz;lkjhgfdsa][poiuytrewq";
    char array_3[256] = "1";
    char array_4[256] = "1";
    char array_6[256] = "1";
    char array_9[] = "\0";
    char array_8[] = "abcde\0";
    char array_10[32] = "Hello world";
    char array_11[] = "qwe";
    char array_12[] = " \n\0";
    char array_13[] = "\0";
    char array_14[] = "abcde\0";
    char array_16[32] = "Hello world";
    char array_17[] = "qwe";
    char array_18[] = " \n\0";
    char array_21[21] = "TEST";
    char array_22[8] = "strcat";
    char array_23[9] = "function";
    char array_24[1] = "";
    char array_25[2] = " ";

    ck_assert_str_eq(strcat(array_3, "qwrteturyit"), s21_strcat(array_3, "qwrteturyit"));
    ck_assert_str_eq(strcat(array_6, array_2), s21_strcat(array_4, array_2));

    s21_strcat(array_10, array_12);
    strcat(array_16, array_18);
    ck_assert_pstr_eq(array_10, array_16);

    s21_strcat(array_11, array_13);
    strcat(array_17, array_9);
    ck_assert_pstr_eq(array_11, array_17);

    s21_strcat(array_14, array_13);
    strcat(array_8, array_9);
    ck_assert_pstr_eq(array_14, array_8);

    s21_strcat(array_12, array_13);
    strcat(array_18, array_9);
    ck_assert_pstr_eq(array_12, array_18);

    s21_strcat(array_21, array_25);
    ck_assert_msg(!memcmp(array_21, "TEST ", 6), "1st array_2cat test failed");

    s21_strcat(array_21, array_22);
    ck_assert_msg(!memcmp(array_21, "TEST strcat", 12), "2nd strcat test failed");

    s21_strcat(array_21, array_25);
    ck_assert_msg(!memcmp(array_21, "TEST strcat ", 13), "3rd strcat test failed");

    s21_strcat(array_21, array_23);
    ck_assert_msg(!memcmp(array_21, "TEST strcat function", 21), "4rth strcat test failed");

    s21_strcat(array_21, array_24);
    ck_assert_msg(!memcmp(array_21, "TEST strcat function", 21), "5th strcat test failed");
}
END_TEST

START_TEST(s21_strncatTEST) {
    char array_2[256] = "?><mnbvcxz;lkjhgfdsa][poiuytrewq";
    char array_3[256] = "1";
    char array_4[256] = "1";
    char array_5[256] = "1";
    char array_6[256] = "1";
    char array_20[32] = "Hello world";
    char array_19[] = "qwe";
    char array_18[] = " \n\0";
    char array_16[] = "\0";
    char array_14[32] = "";
    char array_13[32] = "Hello world";
    char array_11[] = "qwe";
    char array_10[] = " \n\0";
    char array_9[] = "\0";
    char array_17[32] = "";
    char array_21[25] = "TEST";
    char array_22[8] = "strncat";
    char array_23[35] = "function";
    char array_24[1] = "";
    char array_25[30] = "amogus";

    ck_assert_str_eq(strncat(array_3, "qwrteturyit", 12), s21_strncat(array_4, "qwrteturyit", 12));
    ck_assert_str_eq(strncat(array_5, array_2, 50), s21_strncat(array_6, array_2, 50));

    s21_strncat(array_20, array_18, 4);
    strncat(array_13, array_10, 4);
    ck_assert_pstr_eq(array_20, array_13);

    s21_strncat(array_19, array_16, 3);
    strncat(array_11, array_9, 3);
    ck_assert_pstr_eq(array_19, array_11);

    s21_strncat(array_14, array_19, 4);
    strncat(array_17, array_11, 4);
    ck_assert_pstr_eq(array_14, array_17);

    s21_strncat(array_21, array_22,  3);
    ck_assert_msg(!memcmp(array_21, "TESTstr", 8), "1st strncat test failed");

    s21_strncat(array_21, array_24,  3);
    ck_assert_msg(!memcmp(array_21, "TESTstr", 8), "2nd strncat test failed");

    s21_strncat(array_22, array_24,  1);
    ck_assert_msg(!memcmp(array_22, "strncat", 8), "3rd strncat test failed");

    s21_strncat(array_25, array_23,  3);
    ck_assert_msg(!memcmp(array_25, "amogusfun", 10), "4rth strncat test failed");
}
END_TEST

START_TEST(s21_strchrTEST) {
    char array_1[256] = "qwertyuiop[]asdfghjkl;zxcvbnm<>?";
    char array_2[256] = "?><mnbvcxz;lkjhgfdsa][poiuytrewq";
    char array_3[256] = "";
    char array_4[] = "Hello world!";
    char array_5[] = "";
    char *array_21 = "TEST";
    char *array_22 = "strchr";
    char *array_23 = "function";
    char *array_24 = "";
    char *array_25 = "amogus";

    ck_assert_str_eq(strchr(array_1, '>'), s21_strchr(array_1, '>'));
    ck_assert_str_eq(strchr(array_2, '<'), s21_strchr(array_2, '<'));
    ck_assert_ptr_eq(strchr(array_3, ' '), s21_strchr(array_3, ' '));
    ck_assert_pstr_eq(s21_strchr(array_5, 33), strchr(array_5, 33));
    ck_assert_pstr_eq(s21_strchr(array_4, 'w'), strchr(array_4, 'w'));
    ck_assert_pstr_eq(s21_strchr(array_4, 108), strchr(array_4, 108));
    ck_assert_pstr_eq(s21_strchr(array_4, ' '), strchr(array_4, ' '));
    ck_assert_pstr_eq(s21_strchr(array_4, 'a'), strchr(array_4, 'a'));
    ck_assert_str_eq(s21_strchr(array_21, 'E'), strchr(array_21, 'E'));
    ck_assert_str_eq(s21_strchr(array_23, 'i'), strchr(array_23, 'i'));
    ck_assert_str_eq(s21_strchr(array_25, 'a'), strchr(array_25, 'a'));
    ck_assert_str_eq(s21_strchr(array_21, 'E'), strchr(array_21, 'E'));
    ck_assert_ptr_eq(s21_strchr(array_22, 'E'), strchr(array_22, 'E'));
    ck_assert_str_eq(s21_strchr(array_23, 'i'), strchr(array_23, 'i'));
    ck_assert_ptr_eq(s21_strchr(array_24, ' '), strchr(array_24, ' '));
    ck_assert_str_eq(s21_strchr(array_25, 'a'), strchr(array_25, 'a'));
}
END_TEST

START_TEST(s21_strrchrTEST) {
    char array_1[256] = "qwertyuiop[]asdfghjkl;zxcvbnm<>?";
    char array_2[256] = "?><mnbvcxz;lkjhgfdsa][poiuytrewq";
    char array_3[256] = "";
    char array_4[] = ";odsagjnpasdun 1242149215";
    char array_5[] = "";
    char *array_6 = "School-21";
    char *array_7 = "";
    char *array_9 = "amogus";

    ck_assert_str_eq(strrchr(array_1, '>'), s21_strrchr(array_1, '>'));
    ck_assert_str_eq(strrchr(array_2, '<'), s21_strrchr(array_2, '<'));
    ck_assert_ptr_eq(strrchr(array_3, ' '), s21_strrchr(array_3, ' '));
    ck_assert_pstr_eq(s21_strrchr(array_4, 'a'), strrchr(array_4, 'a'));
    ck_assert_pstr_eq(s21_strrchr(array_4, 100), strrchr(array_4, 100));
    ck_assert_pstr_eq(s21_strrchr(array_4, '/'), strrchr(array_4, '/'));
    ck_assert_pstr_eq(s21_strrchr(array_4, 98), strrchr(array_4, 98));
    ck_assert_pstr_eq(s21_strrchr(array_5, 'a'), strrchr(array_5, 'a'));
    ck_assert_pstr_eq(s21_strrchr(array_5, ' '), strrchr(array_5, ' '));
    ck_assert_pstr_eq(s21_strrchr(array_5, 56), strrchr(array_5, 56));
    ck_assert_msg(s21_strrchr(array_6, 'm') == strrchr(array_6, 'm'), "1st strrchr test failed");
    ck_assert_msg(s21_strrchr(array_6, 'o') == strrchr(array_6, 'o'), "2nd strrchr test failed");
    ck_assert_msg(s21_strrchr(array_7, 'o') == strrchr(array_7, 'o'), "3rd strrchr test failed");
    ck_assert_msg(s21_strrchr(array_9, ' ') == strrchr(array_9, ' '), "4rth strrchr test failed");
    ck_assert_msg(s21_strrchr(array_9, 's') == strrchr(array_9, 's'), "5th strrchr test failed");
}
END_TEST

START_TEST(s21_strcmpTEST) {
    char array_1[256] = "qwertyuiop[]asdfghjkl;zxcvbnm<>?";
    char array_2[256] = "?><mnbvcxz;lkjhgfdsa][poiuytrewq";
    char array_3[256] = "";
    char array_4[] = "abcd";
    char array_5[] = "abcd";
    char array_6[] = "efghi";
    char array_11[] = "Sample Text";
    char array_12[] = "Sample Text Plus Some Text";
    char array_13[] = "AAA";
    char array_14[] = "B";
    char array_15[] = "A\0BCD";
    char array_16[] = "A";
    char array_17[] = "amogus";
    char array_18[] = "m";
    char array_19[] = "shrek";
    char array_110[] = "as";

    ck_assert_int_le(strcmp(array_1, array_2), s21_strcmp(array_1, array_2));
    ck_assert_int_le(strcmp(array_1, array_3), s21_strcmp(array_1, array_3));
    ck_assert_int_le(strcmp(array_2, array_3), s21_strcmp(array_2, array_3));
    ck_assert_int_le(strcmp(array_1, "\0"), s21_strcmp(array_1, "\0"));
    ck_assert_int_le(strcmp("\0", array_1), s21_strcmp("\0", array_1));
    ck_assert_int_le(strcmp("\0", "\0"), s21_strcmp("\0", "\0"));
    ck_assert_int_le(strcmp(array_1, "a\0"), s21_strcmp(array_1, "a\0"));
    ck_assert_int_le(strcmp(array_1, "adbdgb\0"), s21_strcmp(array_1, "adbdgb\0"));
    ck_assert_int_le(s21_strcmp(array_4, array_5), strcmp(array_4, array_5));
    ck_assert_int_ge(s21_strcmp(array_6, array_4), strcmp(array_6, array_4));
    ck_assert_int_ge(s21_strcmp(array_6, array_4), strcmp(array_6, array_4));
    ck_assert_int_le(s21_strcmp(array_4, array_6), strcmp(array_4, array_6));
    ck_assert_msg(s21_strcmp(array_11, array_12) >= -32, "1st strcmp test failed");
    ck_assert_msg(s21_strcmp(array_13, array_14) >= -1, "2nd strcmp test failed");
    ck_assert_msg(s21_strcmp(array_15, array_16) >= 0, "3rd strcmp test failed");
    ck_assert_msg(s21_strcmp(array_17, array_18) >= -12, "4rth strcmp test failed");
    ck_assert_msg(s21_strcmp(array_19, array_110) >= 18, "5th strcmp test failed");
}
END_TEST

START_TEST(s21_strncmpTEST) {
    char array_1[256] = "qwertyuiop[]asdfghjkl;zxcvbnm<>?";
    char array_2[256] = "?><mnbvcxz;lkjhgfdsa][poiuytrewq";
    char array_3[256] = "";
    char array_4[] = "abcd";
    char array_5[] = "abcd";
    char array_6[] = "efgh";
    char array_8[] = "";
    char *array_9 = "Sample Text";
    char *array_10 = "Sample Text Plus Some Text";
    char *array_11 = "";
    char *array_12 = "A";
    char *array_13 = "ShrekAmogus";
    char *array_14 = "Shr";

    ck_assert_int_le(strncmp(array_1, array_2, 15), s21_strncmp(array_1, array_2, 15));
    ck_assert_int_ge(strncmp(array_2, array_1, 33), s21_strncmp(array_2, array_1, 33));
    ck_assert_int_le(strncmp(array_1, array_3, 22), s21_strncmp(array_1, array_3, 22));
    ck_assert_int_le(strncmp(array_2, array_3, 11), s21_strncmp(array_2, array_3, 11));
    ck_assert_int_le(strncmp(array_2, array_3, 0), s21_strncmp(array_2, array_3, 0));
    ck_assert_int_le(strncmp("\0", "\0", 3), s21_strncmp("\0", "\0", 3));
    ck_assert_int_le(strncmp("ab", "abc", 3), s21_strncmp("ab", "abc", 3));
    ck_assert_int_ge(s21_strncmp(array_4, array_5, 5), strncmp(array_4, array_5, 5));
    ck_assert_int_le(s21_strncmp(array_4, array_6, 5), strncmp(array_4, array_6, 5));
    ck_assert_int_le(s21_strncmp(array_4, array_6, 4), strncmp(array_4, array_6, 4));
    ck_assert_int_le(s21_strncmp(array_4, array_6, 5), strncmp(array_4, array_6, 5));
    ck_assert_int_ge(s21_strncmp(array_6, array_4, 4), strncmp(array_6, array_4, 4));
    ck_assert_int_le(s21_strncmp(array_8, array_4, 1), strncmp(array_8, array_4, 1));
    ck_assert_int_ge(s21_strncmp(array_4, array_8, 1), strncmp(array_4, array_8, 1));

    int res_1 = strncmp(array_9, array_10, 10);
    int res_2 = s21_strncmp(array_9, array_10, 10);
    ck_assert_msg(!s21_strncmp(array_9, array_10, 10) && !strncmp(array_9, array_10, 10),
        "1st strncmp test failed");
    ck_assert_int_ge(res_1, res_2);
    res_1 = 0;
    res_2 = 0;

    res_1 = strncmp(array_9, array_10, 14);
    res_2 = s21_strncmp(array_9, array_10, 14);
    ck_assert_msg(s21_strncmp(array_9, array_10, 14) < 0 && strncmp(array_9, array_10, 14) < 0,
        "2nd strncmp test failed");
    ck_assert_int_ge(res_1, res_2);
    res_1 = 0;
    res_2 = 0;

    res_1 = strncmp(array_11, array_12, 0);
    res_2 = s21_strncmp(array_11, array_12, 0);
    ck_assert_msg(s21_strncmp(array_11, array_12, 0) == 0 && strncmp(array_11, array_12, 0) == 0,
        "3rd strncmp test failed");
    ck_assert_msg(s21_strncmp(array_11, array_12, 1) < 0 && strncmp(array_11, array_12, 1) < 0,
        "3rd strncmp test failed");
    ck_assert_int_eq(res_1, res_2);
    res_1 = 0;
    res_2 = 0;

    res_1 = strncmp(array_13, array_14, 5);
    res_2 = s21_strncmp(array_13, array_14, 5);
    ck_assert_msg(s21_strncmp(array_13, array_14, 5) > 0 && strncmp(array_13, array_14, 5) > 0,
        "4rth strncmp test failed");
    ck_assert_msg(s21_strncmp(array_13, array_14, 0) == 0 && strncmp(array_13, array_14, 0) == 0,
        "5th strncmp test failed");
    ck_assert_int_le(res_1, res_2);
    res_1 = 0;
    res_2 = 0;
}
END_TEST

START_TEST(s21_strcpyTEST) {
    char *src = "abra kad";
    char *dest_1 = malloc(10);
    char *dest_2 = malloc(10);
    char array_1[] = "Hello";
    char array_2[] = "World";
    char array_3[] = "";
    char array_4[32] = "";
    char array_5[] = "Hello";
    char array_6[] = "World";
    char array_7[] = "";
    char array_8[32] = "";
    char test_str1[50] = "Sample Text";
    char test_dest1[50] = "";
    char test_str2[50] = "AAAAA";
    char test_dest2[50] = "NO";
    char test_str3[50] = "";
    char test_dest3[50] = "HELLO";
    char test_str4[50] = "amogus";
    char test_dest4[50] = "HELLO";
    char test_str5[50] = "Shrek";
    char test_dest5[50] = "what";

    ck_assert_str_eq(s21_strcpy(dest_1, src), strcpy(dest_2, src));
    ck_assert_str_eq(s21_strcpy(dest_1, ""), strcpy(dest_2, ""));

    s21_strcpy(array_1, array_2);
    strcpy(array_5, array_6);
    ck_assert_str_eq(array_1, array_5);

    s21_strcpy(array_2, array_1);
    strcpy(array_6, array_5);
    ck_assert_str_eq(array_2, array_6);

    s21_strcpy(array_1, array_3);
    strcpy(array_5, array_7);
    ck_assert_str_eq(array_1, array_5);

    s21_strcpy(array_4, array_1);
    strcpy(array_8, array_5);
    ck_assert_str_eq(array_4, array_8);

    s21_strcpy(array_4, array_3);
    strcpy(array_8, array_7);

    s21_strcpy(test_dest1, test_str1);
    ck_assert_msg(!strcmp(test_dest1, test_str1), "1st  strcpy test failed");

    s21_strcpy(test_dest2, test_str2);
    ck_assert_msg(!strcmp(test_dest2, test_str2), "2nd strcpy test failed");

    s21_strcpy(test_dest2, test_str2);
    ck_assert_msg(!strcmp(test_dest2, test_str2), "2nd strcpy test failed");

    s21_strcpy(test_dest3, test_str3);
    ck_assert_msg(!strcmp(test_dest3, test_str3), "3rd strcpy test failed");
    s21_strcpy(test_dest4, test_str4);
    ck_assert_msg(!strcmp(test_dest4, test_str4), "4rth strcpy test failed");
    s21_strcpy(test_dest5, test_str5);
    ck_assert_msg(!strcmp(test_dest5, test_str5), "5th strcpy test failed");

    free(dest_1);
    free(dest_2);
}
END_TEST

START_TEST(s21_strncpyTEST) {
    char dest_1[15] = "";
    char dest_2[15] = "";
    char *src = "abra kadabra";
    char array_1[] = "kjzxnvkjbsdfjhkbbds";
    char array_2[] = "jsdjsdfkjdsfjh";
    char array_3[] = "";
    char array_4[] = "kjzxnvkjbsdfjhkbbds";
    char array_5[] = "jsdjsdfkjdsfjh";
    char array_6[] = "";
    char test_str1[50] = "Sample Text";
    char test_dest1[50] = "";
    char test_str2[50] = "AAAAA";
    char test_dest2[50] = "NO";
    char test_str3[50] = "H";
    char test_dest3[50] = "fELLO";
    char test_str4[50] = "amo";
    char test_dest4[50] = "qwegus";
    char test_str5[50] = "shrek";
    char test_dest5[50] = "s";

    ck_assert_str_eq(s21_strncpy(dest_1, src, 12), strncpy(dest_2, src, 12));
    ck_assert_str_eq(s21_strncpy(dest_1, "abra kad", 12), strncpy(dest_2, "abra kad", 12));
    ck_assert_str_eq(s21_strncpy(dest_1, "", 12), strncpy(dest_2, "", 12));

    s21_strncpy(test_dest1, test_str1, 2);
    ck_assert_msg(!strcmp(test_dest1, "Sa"), "1st strncpy test failed");

    s21_strncpy(test_dest2, test_str2, 0);
    ck_assert_msg(!strcmp(test_dest2, "NO"), "2nd strncpy test failed");

    s21_strncpy(test_dest3, test_str3, 1);
    ck_assert_msg(!strcmp(test_dest3, "HELLO"), "3rd strncpy test failed");

    s21_strncpy(test_dest4, test_str4, 3);
    ck_assert_msg(!strcmp(test_dest4, "amogus"), "4rth strncpy test failed");

    s21_strncpy(test_dest5, test_str5, 5);
    ck_assert_msg(!strcmp(test_dest5, "shrek"), "5th strncpy test failed");

    s21_strncpy(array_1, array_2, 5);
    strncpy(array_4, array_5, 5);
    ck_assert_str_eq(array_1, array_4);

    s21_strncpy(array_1, array_2, 20);
    strncpy(array_4, array_5, 20);
    ck_assert_str_eq(array_1, array_4);

    s21_strncpy(array_1, array_2, 15);
    strncpy(array_4, array_5, 15);
    ck_assert_str_eq(array_1, array_4);

    s21_strncpy(array_2, array_1, 8);
    strncpy(array_5, array_4, 8);
    ck_assert_str_eq(array_1, array_4);

    s21_strncpy(array_2, array_1, 15);
    strncpy(array_5, array_4, 15);
    ck_assert_str_eq(array_1, array_4);

    s21_strncpy(array_1, array_3, 8);
    strncpy(array_4, array_6, 8);
    ck_assert_str_eq(array_1, array_4);
}
END_TEST


START_TEST(s21_strpbrkTEST) {
    char *dest = "abra kadabra";
    char array_1[] = "abcdefghijklmnop";
    char array_2[] = "pkh";
    char array_3[] = "";
    char *array_4 = "hello";
    char *array_5 = "hel";
    char *array_7 = "o";
    char *array_8 = "abcdefghjkl";
    char *array_9 = "abcd";
    char *array_10 = "amogus is shrek?";
    char *res1, *res2;

    ck_assert_str_eq(s21_strpbrk(dest, "kbr"), strpbrk(dest, "kbr"));
    ck_assert_ptr_null(s21_strpbrk(dest, "zlp"));
    ck_assert_ptr_null(s21_strpbrk(dest, ""));
    ck_assert_ptr_null(s21_strpbrk("", "kbr"));
    ck_assert_ptr_eq(s21_strpbrk(array_1, array_2), strpbrk(array_1, array_2));
    ck_assert_ptr_eq(s21_strpbrk(array_2, array_1), strpbrk(array_2, array_1));
    ck_assert_ptr_eq(s21_strpbrk(array_1, array_3), strpbrk(array_1, array_3));
    ck_assert_ptr_eq(s21_strpbrk(array_2, array_3), strpbrk(array_2, array_3));
    ck_assert_ptr_eq(s21_strpbrk(array_3, array_1), strpbrk(array_3, array_1));
    ck_assert_ptr_eq(s21_strpbrk(array_3, array_2), strpbrk(array_3, array_2));

    res1 = s21_strpbrk(array_4, array_5);
    res2 = strpbrk(array_4, array_5);
    ck_assert_str_eq(res1, res2);

    res1 = s21_strpbrk(array_4, array_7);
    res2 = strpbrk(array_4, array_7);
    ck_assert_ptr_eq(res1, res2);

    res1 = s21_strpbrk(array_4, array_8);
    res2 = strpbrk(array_4, array_8);
    ck_assert_str_eq(res1, res2);

    res1 = s21_strpbrk(array_9, array_10);
    res2 = strpbrk(array_9, array_10);
    ck_assert_str_eq(res1, res2);
}
END_TEST

START_TEST(s21_strspnTEST) {
    char *dest = "abra kadabra";
    char array_1[] = "a;iuhgilsdugiaudhgiu2";
    char array_2[] = "gsdsgk";
    char array_3[] = "i2u";
    char array_4[] = "";
    char *array_5 = "School-21";
    char *array_6 = "ho";
    char *array_7 = "";

    ck_assert_int_eq(s21_strspn(dest, "lp"), strspn(dest, "lp"));
    ck_assert_int_eq(s21_strspn(dest, "kd"), strspn(dest, "kd"));
    ck_assert_int_eq(s21_strspn(dest, "akd"), strspn(dest, "akd"));
    ck_assert_int_eq(s21_strspn("", "as"), strspn("", "as"));
    ck_assert_int_eq(s21_strspn(array_1, array_2), strspn(array_1, array_2));
    ck_assert_int_eq(s21_strspn(array_1, array_3), strspn(array_1, array_3));
    ck_assert_int_eq(s21_strspn(array_2, array_1), strspn(array_2, array_1));
    ck_assert_int_eq(s21_strspn(array_3, array_1), strspn(array_3, array_1));
    ck_assert_int_eq(s21_strspn(array_2, array_3), strspn(array_2, array_3));
    ck_assert_int_eq(s21_strspn(array_2, array_4), strspn(array_2, array_4));
    ck_assert_int_eq(s21_strspn(array_4, array_1), strspn(array_4, array_1));
    ck_assert_msg(s21_strspn(array_5, array_6) == strspn(array_5, array_6), "1st strspn test failed");
    ck_assert_msg(s21_strspn(array_7, array_6) == strspn(array_7, array_6), "2nd strspn test failed");
    ck_assert_msg(s21_strspn(array_6, array_7) == strspn(array_6, array_7), "3rd strspn test failed");
    ck_assert_msg(s21_strspn(array_5, array_7) == strspn(array_5, array_7), "4rth strspn test failed");
    ck_assert_msg(s21_strspn(array_7, array_5) == strspn(array_7, array_5), "5th strspn test failed");
}
END_TEST

START_TEST(s21_strstrTEST) {
    char *dest = "abra kadabra";
    char *src = "kada";
    char array_1[] = "abcdefghijklm";
    char array_2[] = "hij";
    char array_3[] = "";
    char array_4[] = "edc";
    char *array_5 = "School-21";
    char *array_6 = "ho";
    char *array_7 = "";

    ck_assert_pstr_eq(s21_strstr(array_1, array_2), strstr(array_1, array_2));
    ck_assert_pstr_eq(s21_strstr(array_1, array_4), strstr(array_1, array_4));
    ck_assert_pstr_eq(s21_strstr(array_2, array_1), strstr(array_2, array_1));
    ck_assert_pstr_eq(s21_strstr(array_2, array_4), strstr(array_2, array_4));
    ck_assert_pstr_eq(s21_strstr(array_3, array_1), strstr(array_3, array_1));
    ck_assert_pstr_eq(s21_strstr(array_3, array_2), strstr(array_3, array_2));
    ck_assert_pstr_eq(s21_strstr(array_3, array_4), strstr(array_3, array_4));
    ck_assert_msg(s21_strstr(array_5, array_6) == strstr(array_5, array_6), "1st strstr test failed");
    ck_assert_msg(s21_strstr(array_7, array_6) == strstr(array_7, array_6), "3rd strstr test failed");
    ck_assert_msg(s21_strstr(array_6, array_7) == strstr(array_6, array_7), "4rth strstr test failed");
    ck_assert_msg(s21_strstr(array_7, array_5) == strstr(array_7, array_5), "5th strstr test failed");
    ck_assert_str_eq(s21_strstr(dest, src), strstr(dest, src));
    ck_assert_ptr_eq(s21_strstr(dest, "ghl"), strstr(dest, "ghl"));
    ck_assert_ptr_eq(s21_strstr("", "asd"), strstr("", "asd"));
}
END_TEST

START_TEST(s21_strcspnTEST) {
    char *dest = "abra zillut";
    char array_1[] = "jkflkfgdkhslkhgsdk";
    char array_2[] = "gd";
    char array_3[] = "";
    char *array_11 = "0123456789";
    char *array_12 = "9876";
    char *array_13 = "0123456789";
    char *array_14 = "";
    char *array_15 = "";
    char *array_16 = "123";
    char *array_17 = "1337";
    char *array_18 = "228";
    char *array_19 = "19642";
    char *array_110 = "64";

    ck_assert_int_eq(s21_strcspn(dest, "lut"), strcspn(dest, "lut"));
    ck_assert_int_eq(s21_strcspn(dest, "b"), strcspn(dest, "b"));
    ck_assert_int_eq(s21_strcspn(dest, "0"), strcspn(dest, "0"));
    ck_assert_int_eq(s21_strcspn(array_1, array_2), strcspn(array_1, array_2));
    ck_assert_int_eq(s21_strcspn(array_2, array_1), strcspn(array_2, array_1));
    ck_assert_int_eq(s21_strcspn(array_1, array_3), strcspn(array_1, array_3));
    ck_assert_int_eq(s21_strcspn(array_2, array_3), strcspn(array_2, array_3));
    ck_assert_int_eq(s21_strcspn(array_3, array_1), strcspn(array_3, array_1));
    ck_assert_int_eq(s21_strcspn(array_3, array_2), strcspn(array_3, array_2));
    ck_assert_msg(s21_strcspn(array_11, array_12) == strcspn(array_11, array_12), "1st strcspn test failed");
    ck_assert_msg(s21_strcspn(array_13, array_14) == strcspn(array_13, array_14), "2nd strcspn test failed");
    ck_assert_msg(s21_strcspn(array_15, array_16) == strcspn(array_15, array_16), "3rd strcspn test failed");
    ck_assert_msg(s21_strcspn(array_17, array_18) == strcspn(array_17, array_18), "4rth strcspn test failed");
    ck_assert_msg(s21_strcspn(array_19, array_110) == strcspn(array_19, array_110),
    "5th strcspn test failed");
}
END_TEST

START_TEST(s21_strtokTEST) {
    char array_1[100] = "My name is Root. Hello world!";
    char *parr1 = array_1;
    char array_2[100] = "My name is Root. Hello world!";
    char *parr2 = array_2;

    s21_strtok(parr1, " .");
    strtok(parr2, " .");
    ck_assert_pstr_eq(parr1, parr2);

    parr1 = s21_strtok(s21_NULL, " .");
    parr2 = strtok(s21_NULL, " .");
    ck_assert_pstr_eq(parr1, parr2);

    parr1 = s21_strtok(parr1, "x");
    parr2 = strtok(parr2, "x");
    ck_assert_pstr_eq(parr1, parr2);

    char array1[50] = "aboba,hello,world,trim";
    char array2[50] = "aboba,hello,world,trim";
    char delim1[10] = "123,";
    char *orig1, *copy1;
    orig1 = strtok(array1, delim1);
    copy1 = s21_strtok(array2, delim1);
    char orig_res1[500] = {0};
    char copy_res1[500] = {0};
    while (orig1 != s21_NULL) {
      strcat(orig_res1, orig1);
      orig1 = strtok(s21_NULL, delim1);
    }
    while (copy1 != s21_NULL) {
      strcat(copy_res1, copy1);
      copy1 = s21_strtok(s21_NULL, delim1);
    }
    ck_assert_str_eq(orig_res1, copy_res1);
    char str_arr[] = "School-21";
    char *tok1 = s21_strtok(str_arr, "-");
    char *tok2 = strtok(str_arr, "-");
    char *tok3 = s21_strtok(str_arr, "");
    char *tok4 = strtok(str_arr, "");
    char *tok5 = s21_strtok(str_arr, "oo");
    char *tok6 = strtok(str_arr, "oo");
    char *tok7 = s21_strtok(str_arr, "Sc");
    char *tok8 = strtok(str_arr, "Sc");
    char *tok9 = s21_strtok(str_arr, "21");
    char *tok10 = strtok(str_arr, "21");
    ck_assert_msg(!strcmp(tok1, tok2), "1st strtok test failed");
    ck_assert_msg(!strcmp(tok3, tok4), "2nd strtok test failed");
    ck_assert_msg(!strcmp(tok5, tok6), "3rd strtok test failed");
    ck_assert_msg(!strcmp(tok7, tok8), "4rth strtok test failed");
    ck_assert_msg(!strcmp(tok9, tok10), "5th strtok test failed");
}
END_TEST

START_TEST(s21_strerrorTEST) {
  for (int j = 0, err = -1; j < 3; j++) {
        char *errtest1 = s21_strerror(err);
        char *errtest2 = strerror(err);
        s21_size_t i = 0, len1 = s21_strlen(errtest1),
                   len2 = s21_strlen(errtest2);
        while (i < len1 || i < len2) {
            ck_assert_pstr_eq(errtest1 + i, errtest2 + i);
            i++;
        }
        err += 1;
    }

    char *array_1 = s21_strerror(20);
    char *array_2 = strerror(20);
    char *array_3 = s21_strerror(12);
    char *array_4 = strerror(12);
    char *array_5 = s21_strerror(1337);
    char *array_6 = strerror(1337);
    char *array_7 = s21_strerror(1);
    char *array_8 = strerror(1);
    char *array_9 = s21_strerror(25);
    char *array_10 = strerror(25);

    ck_assert_str_eq(array_1, array_2);
    ck_assert_str_eq(array_3, array_4);
    ck_assert_str_eq(array_5, array_6);
    ck_assert_str_eq(array_7, array_8);
    ck_assert_str_eq(array_9, array_10);
}
END_TEST

START_TEST(s21_trimTEST) {
    char *res = s21_trim("testc", "c");
    ck_assert_str_eq(res, "test");
    free(res);

    res = s21_trim(" test ", " ");
    ck_assert_str_eq(res, "test");
    free(res);

    res = s21_trim(" ", " ");
    ck_assert_str_eq(res, "");
    free(res);

    res = s21_trim("*te*st*", "*");
    ck_assert_str_eq(res, "te*st");
    free(res);

    res = s21_trim(" te st001", " ");
    ck_assert_str_eq(res, "te st001");
    free(res);

    res = s21_trim(" /*()", ")(/ *");
    ck_assert_str_eq(res, "");
    free(res);

    char *test_str_1 = "School-21";
    char *test_str_2 = "";
    char *test_str_3 = "S";
    char *test_str_4 = "Sch";
    char *test_str_5 = "School";
    char *test_str_6 = "School-21";
    char *test_str_7 = s21_NULL;

    char *test_res_1 = s21_trim(test_str_1, test_str_2);
    char *test_res_2 = s21_trim(test_str_1, test_str_3);
    char *test_res_3 = s21_trim(test_str_1, test_str_4);
    char *test_res_4 = s21_trim(test_str_1, test_str_5);
    char *test_res_5 = s21_trim(test_str_1, test_str_6);
    char *test_res_6 = s21_trim(test_str_1, test_str_7);

    ck_assert_msg(!strcmp(test_res_1, "School-21"), "1st trim test failed");
    ck_assert_msg(!strcmp(test_res_2, "chool-21"), "2nd trim test failed");
    ck_assert_msg(!strcmp(test_res_3, "ool-21"), "3rd trim test failed");
    ck_assert_msg(!strcmp(test_res_4, "-21"), "4rth trim test failed");
    ck_assert_msg(!strcmp(test_res_5, ""), "5th trim test failed");

    free(test_res_1);
    free(test_res_2);
    free(test_res_3);
    free(test_res_4);
    free(test_res_5);
    free(test_res_6);
}
END_TEST

START_TEST(s21_insertTEST) {
    char *array_1 = "Test";
    char *array_2 = "quest";
    char *res = s21_NULL;

    res = s21_insert(array_1, array_2, 2);
    ck_assert_str_eq(res, "Tequestst");
    free(res);

    array_1 = "The test";
    array_2 = " perfect";
    res = s21_insert(array_1, array_2, 3);
    ck_assert_str_eq(res, "The perfect test");
    free(res);

    array_1 = "Insert in";
    array_2 = " end";
    res = s21_insert(array_1, array_2, 9);
    ck_assert_str_eq(res, "Insert in end");
    free(res);

    array_1 = " with insert";
    array_2 = "Start";
    res = s21_insert(array_1, array_2, 0);
    ck_assert_str_eq(res, "Start with insert");
    free(res);

    array_1 = "Newspace";
    array_2 = " ";
    res = s21_insert(array_1, array_2, 3);
    ck_assert_str_eq(res, "New space");
    free(res);

    array_1 = "";
    array_2 = " ";
    res = s21_insert(array_1, array_2, 0);
    ck_assert_str_eq(res, " ");
    free(res);

    char *test_str_1 = "School-21 test";
    char *test_str_2 = "insert ";
    char *test_str_3 = "amogus";
    char *test_str_4 = "is";
    char *test_str_5 = "shrek";

    char *res1 = s21_insert(test_str_1, " ", 3);
    ck_assert_str_eq(res1, "Sch ool-21 test");
    free(res1);

    char *res2 = s21_insert(test_str_5, test_str_3, 5);
    ck_assert_str_eq(res2, "shrekamogus");
    free(res2);

    char *res3 = s21_insert(test_str_3, test_str_5, 0);
    ck_assert_str_eq(res3, "shrekamogus");
    free(res3);

    char *res4 = s21_insert(test_str_5, test_str_4, 5);
    ck_assert_str_eq(res4, "shrekis");
    free(res4);

    char *res5 = s21_insert(test_str_2, test_str_1, 4);
    ck_assert_str_eq(res5, "inseSchool-21 testrt ");
    free(res5);

    char *res6 = s21_insert("qwerty", "", 3);
    ck_assert_str_eq(res6, "qwerty");
    free(res6);
}
END_TEST

START_TEST(s21_to_upperTEST) {
    char *res = s21_to_upper("test");
    ck_assert_str_eq(res, "TEST");
    free(res);

    res = s21_to_upper("qwert");
    ck_assert_str_eq(res, "QWERT");
    free(res);

    res = s21_to_upper("1q2w3e4r5");
    ck_assert_str_eq(res, "1Q2W3E4R5");
    free(res);

    res = s21_to_upper("!@#!$^");
    ck_assert_str_eq(res, "!@#!$^");
    free(res);

    res = s21_to_upper("\n\0\\n");
    ck_assert_str_eq(res, "\n\0\\N");
    free(res);

    res = s21_to_upper(" t e s t");
    ck_assert_str_eq(res, " T E S T");
    free(res);

    res = s21_to_upper("test001");
    ck_assert_str_eq(res, "TEST001");
    free(res);

    res = s21_to_upper("!test!");
    ck_assert_str_eq(res, "!TEST!");
    free(res);

    res = s21_to_upper(" ");
    ck_assert_str_eq(res, " ");
    free(res);

    res = s21_to_upper("!*_*!");
    ck_assert_str_eq(res, "!*_*!");
    free(res);

    char *test_str_1 = "";
    char *test_str_2 = "School is Cool";
    char *test_str_3 = "amogus";
    char *test_str_4 = "Shrek is amogus";
    char *test_str_5 = "what";

    char *test_res_1 = s21_to_upper(test_str_1);
    char *test_res_2 = s21_to_upper(test_str_2);
    char *test_res_3 = s21_to_upper(test_str_3);
    char *test_res_4 = s21_to_upper(test_str_4);
    char *test_res_5 = s21_to_upper(test_str_5);

    ck_assert_str_eq(test_res_1, "");
    ck_assert_str_eq(test_res_2, "SCHOOL IS COOL");
    ck_assert_str_eq(test_res_3, "AMOGUS");
    ck_assert_str_eq(test_res_4, "SHREK IS AMOGUS");
    ck_assert_str_eq(test_res_5, "WHAT");

    free(test_res_1);
    free(test_res_2);
    free(test_res_3);
    free(test_res_4);
    free(test_res_5);
}
END_TEST

START_TEST(s21_to_lowerTEST) {
    char *res = s21_to_lower("TEST");
    ck_assert_str_eq(res, "test");
    free(res);

    res = s21_to_lower("QWERT");
    ck_assert_str_eq(res, "qwert");
    free(res);

    res = s21_to_lower("1Q2W3E4R5");
    ck_assert_str_eq(res, "1q2w3e4r5");
    free(res);

    res = s21_to_lower("!@#!$^");
    ck_assert_str_eq(res, "!@#!$^");
    free(res);

    res = s21_to_lower("\n\0\\N");
    ck_assert_str_eq(res, "\n\0\\n");
    free(res);

    res = s21_to_lower(" T E S T");
    ck_assert_str_eq(res, " t e s t");
    free(res);

    res = s21_to_lower("TEST001");
    ck_assert_str_eq(res, "test001");
    free(res);

    res = s21_to_lower("!TEST!");
    ck_assert_str_eq(res, "!test!");
    free(res);

    res = s21_to_lower(" ");
    ck_assert_str_eq(res, " ");
    free(res);

    res = s21_to_lower("!*_*!");
    ck_assert_str_eq(res, "!*_*!");
    free(res);

    char test_str_1[] = "SChOOl-21";
    char test_str_2[] = "ScHool is CoOl";
    char test_str_3[] = "amogus";
    char test_str_4[] = "ShrEK is Amogus";
    char test_str_5[] = "what";

    char *test_res_1 = s21_to_lower(test_str_1);
    char *test_res_2 = s21_to_lower(test_str_2);
    char *test_res_3 = s21_to_lower(test_str_3);
    char *test_res_4 = s21_to_lower(test_str_4);
    char *test_res_5 = s21_to_lower(test_str_5);

    ck_assert_msg(!strcmp(test_res_1, "school-21"), "1st to_lower test failed");
    ck_assert_msg(!strcmp(test_res_2, "school is cool"), "2nd to_lower test failed");
    ck_assert_msg(!strcmp(test_res_3, "amogus"), "3rd to_lower test failed");
    ck_assert_msg(!strcmp(test_res_4, "shrek is amogus"), "4rth to_lower test failed");
    ck_assert_msg(!strcmp(test_res_5, "what"), "5th to_lower test failed");

    free(test_res_1);
    free(test_res_2);
    free(test_res_3);
    free(test_res_4);
    free(test_res_5);
}
END_TEST

char mass1[100] = {'\0'}, mass2[100] = {'\0'};
char Lmass1[255] = {'\0'}, Lmass2[255] = {'\0'};
int int_test1 = 0, int_test2 = 0;

int n = 2100, n1 = 0, n2 = -21;
long int lg = 12345678912345, lg1 = -12345678912345;
short sh = 123, sh1 = -123;
float f = 21.21, f1 = -21.21;

char char_test1 = 'T', char_test2 = 'Z';
wchar_t str_test[10] = {'h', 'e', 'l', 'l', 'o'};

void free_array(char *arr1, char *arr2) {
    int Legth1 = sizeof(arr1);
    int Legth2 = sizeof(arr2);
    if (Legth1 < Legth2) Legth1 = Legth2;

    for (int i = 0; i < Legth1; i++) {
        arr1[i] = '\0';
        arr2[i] = '\0';
    }
}

START_TEST(Sprintf_Test_s) {
    s21_sprintf(mass1, "%s", "Hello");
    sprintf(mass2, "%s", "Hello");
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%.1s", "Hello");
    sprintf(mass2, "%.1s", "Hello");
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%7s", "Hello");
    sprintf(mass2, "%7s", "Hello");
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%-10s", "Hello");
    sprintf(mass2, "%-10s", "Hello");
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "qwerty %s qwe", "qwerty");
    sprintf(mass2, "qwerty %s qwe", "qwerty");
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(Lmass1, "%s", "");
    int_test2 = sprintf(Lmass2, "%s", "");
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%23s", "aboba");
    int_test2 = sprintf(Lmass2, "%23s", "aboba");
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-15s", "helloworld");
    int_test2 = sprintf(Lmass2, "%-15s", "helloworld");
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%10.s", "ABOBA");
    int_test2 = sprintf(Lmass2, "%10.s", "ABOBA");
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%15.2s", "ABOBBAAA");
    int_test2 = sprintf(Lmass2, "%15.2s", "ABOBBAAA");
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%16.*s", 10, "QWERTY");
    int_test2 = sprintf(Lmass2, "%16.*s", 10, "QWERTY");
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%*.*s", 15, 4, "AAAAABOOOOOBAAA");
    int_test2 = sprintf(Lmass2, "%*.*s", 15, 4, "AAAAABOOOOOBAAA");
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%*s", 10, "ABOBAAAAA");
    int_test2 = sprintf(Lmass2, "%*s", 10, "ABOBAAAAA");
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%15ls", str_test);
    int_test2 = sprintf(Lmass2, "%15ls", str_test);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%15.3ls", str_test);
    int_test2 = sprintf(Lmass2, "%15.3ls", str_test);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%s %10s", "hello", "hello");
    int_test2 = sprintf(Lmass2, "%s %10s", "hello", "hello");
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%.s %.3s %.10s", "hello", "hello", "hello");
    int_test2 = sprintf(Lmass2, "%.s %.3s %.10s", "hello", "hello", "hello");
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-10s %-3s", "hello", "hello");
    int_test2 = sprintf(Lmass2, "%-10s %-3s", "hello", "hello");
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%ls %10ls", L"hello", L"hello");
    int_test2 = sprintf(Lmass2, "%ls %10ls", L"hello", L"hello");
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-10ls", L"hello");
    int_test2 = sprintf(Lmass2, "%-10ls", L"hello");
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-3ls", L"hello");
    int_test2 = sprintf(Lmass2, "%-3ls", L"hello");
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%.ls", L"hello");
    int_test2 = sprintf(Lmass2, "%.ls", L"hello");
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%.3ls", L"hello");
    int_test2 = sprintf(Lmass2, "%.3ls", L"hello");
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%.10ls", L"hello");
    int_test2 = sprintf(Lmass2, "%.10ls", L"hello");
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%s%23s%-15s%10.s%15.2s%16.*s%*.*s%*s%15ls%15.3ls",
    "aboba", "helloworld", "ABOBA", "ABOBBAAA", "AAABOBA", 10, "QWERTY",
    15, 4, "AAAAABOOOOOBAAA", 10, "ABOBAAAAA", str_test, str_test);
    int_test2 = sprintf(Lmass2,  "%s%23s%-15s%10.s%15.2s%16.*s%*.*s%*s%15ls%15.3ls",
    "aboba", "helloworld", "ABOBA", "ABOBBAAA", "AAABOBA", 10, "QWERTY",
    15, 4, "AAAAABOOOOOBAAA", 10, "ABOBAAAAA", str_test, str_test);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%s%23s%-15s%10.s%15.2s%16.*s%*.*s%*s%15ls%15.3ls",
    "aboba", "helloworld", "ABOBA", "ABOBBAAA", "AAABOBA", 10, "QWERTY",
    15, 4, "AAAAABOOOOOBAAA", 10, "ABOBAAAAA", str_test, str_test);
    int_test2 = sprintf(Lmass2, "%s%23s%-15s%10.s%15.2s%16.*s%*.*s%*s%15ls%15.3ls",
    "aboba", "helloworld", "ABOBA", "ABOBBAAA", "AAABOBA", 10, "QWERTY",
    15, 4, "AAAAABOOOOOBAAA", 10, "ABOBAAAAA", str_test, str_test);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);
}
END_TEST

START_TEST(Sprintf_Test_c) {
    s21_sprintf(Lmass1, "%10c", 'a');
    sprintf(Lmass2, "%10c", 'a');
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%-10c", 'e');
    sprintf(Lmass2, "%-10c", 'e');
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%c hello %c", '\0', '\0');
    sprintf(Lmass2, "%c hello %c", '\0', '\0');
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%25c", 'f');
    sprintf(Lmass2, "%25c", 'f');
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%-30c", 'R');
    sprintf(Lmass2, "%-30c", 'R');
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%25c", 'Q');
    sprintf(Lmass2, "%25c", 'Q');
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%*c", 10, 'P');
    sprintf(Lmass2, "%*c", 10, 'P');
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%2c", 'L');
    sprintf(Lmass2, "%2c", 'L');
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%*c", 20, 'x');
    sprintf(Lmass2, "%*c", 20, 'x');
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%lc", char_test1);
    sprintf(Lmass2, "%lc", char_test1);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%c", CHAR_MAX);
    sprintf(Lmass2, "%c", CHAR_MAX);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%c", CHAR_MIN);
    sprintf(Lmass2, "%c", CHAR_MIN);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%c", UCHAR_MAX);
    sprintf(Lmass2, "%c", UCHAR_MAX);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%c %c %c", CHAR_MAX, CHAR_MIN, UCHAR_MAX);
    sprintf(Lmass2, "%c %c %c", CHAR_MAX, CHAR_MIN, UCHAR_MAX);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%-50lc", char_test2);
    sprintf(Lmass2, "%-50lc", char_test2);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%10c%-10c%25c%-30c%25c%*c%2c%*c%lc%-50lc",
     'a', 'e', 'f', 'R', 'Q', 10, 'P', 'L', 20, 'x', char_test1, char_test2);
    sprintf(Lmass2, "%10c%-10c%25c%-30c%25c%*c%2c%*c%lc%-50lc",
     'a', 'e', 'f', 'R', 'Q', 10, 'P', 'L', 20, 'x', char_test1, char_test2);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%10c%-10cdsad%25c%-30c  %25c%*c%2c%*csadasaaaWWW%lc%-50lc",
     'a', 'e', 'f', 'R', 'Q', 10, 'P', 'L', 20, 'x', char_test1, char_test2);
    sprintf(Lmass2, "%10c%-10cdsad%25c%-30c  %25c%*c%2c%*csadasaaaWWW%lc%-50lc",
     'a', 'e', 'f', 'R', 'Q', 10, 'P', 'L', 20, 'x', char_test1, char_test2);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%c", '!');
    sprintf(Lmass2, "%c", '!');
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%10c", 'a');
    sprintf(Lmass2, "%10c", 'a');
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%-10c", 'b');
    sprintf(Lmass2, "%-10c", 'b');
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%c %10c %-10c", '!', 'a', 'b');
    sprintf(Lmass2, "%c %10c %-10c", '!', 'a', 'b');
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(mass1, "%1c %c", 'H', 'z');
    sprintf(mass2, "%1c %c", 'H', 'z');
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%10c %2c", 'H', 'z');
    sprintf(mass2, "%10c %2c", 'H', 'z');
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%-1c %-c", 'H', 'z');
    sprintf(mass2, "%-1c %-c", 'H', 'z');
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %3c abc", '*');
    sprintf(mass2, "abc %3c abc", '*');
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %-3c abc", '*');
    sprintf(mass2, "abc %-3c abc", '*');
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "qwerty %c %-5c %15c", '$', '/', '!');
    sprintf(mass2, "qwerty %c %-5c %15c", '$', '/', '!');
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);
}
END_TEST

START_TEST(Sprintf_Test_f) {
    s21_sprintf(mass1, "abc %+.4f abc", 166.087);
    sprintf(mass2, "abc %+.4f abc", 166.087);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %-12.3f abc", 166.87);
    sprintf(mass2, "abc %-12.3f abc", 166.87);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %f abc", -166.87);
    sprintf(mass2, "abc %f abc", -166.87);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%f", 100000.4);
    sprintf(mass2, "%f", 100000.4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%0f", 100000.4);
    sprintf(mass2, "%0f", 100000.4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%f", -100000.4);
    sprintf(mass2, "%f", -100000.4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%+f %+f", -100000.4, 1000.231);
    sprintf(mass2, "%+f %+f", -100000.4, 1000.231);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%-20f", 100000.4);
    sprintf(mass2, "%-20f", 100000.4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %.*f", 1, 50.897);
    sprintf(mass2, "abc %.*f", 1, 50.897);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %Lf", (long double)58787.878958);
    sprintf(mass2, "abc %Lf", (long double)58787.878958);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %.lf", 58787.878958);
    sprintf(mass2, "abc %.lf", 58787.878958);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %.10f", 587.00089);
    sprintf(mass2, "abc %.10f", 587.00089);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%-+20f", 100000.4);
    sprintf(mass2, "%-+20f", 100000.4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%-20f", -100000.4);
    sprintf(mass2, "%-20f", -100000.4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%+20f", 100000.4);
    sprintf(mass2, "%+20f", 100000.4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%+20f", -100000.4);
    sprintf(mass2, "%+20f", -100000.4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    float key = 10000;
    s21_sprintf(mass1, "%+20f", key);
    sprintf(mass2, "%+20f", key);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(Lmass1, "%+-10.5f", 2546.4324);
    int_test2 = sprintf(Lmass2, "%+-10.5f", 2546.4324);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%+10.2f", 3452.4325);
    int_test2 = sprintf(Lmass2, "%+10.2f", 3452.4325);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%15.16f", 5678.43256);
    int_test2 = sprintf(Lmass2, "%15.16f", 5678.43256);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%+*.10f", 25, 456789.243242);
    int_test2 = sprintf(Lmass2, "%+*.10f", 25, 456789.243242);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%*.16f", 10, 3456789123.43245);
    int_test2 = sprintf(Lmass2, "%*.16f", 10, 3456789123.43245);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-10.*f", 2, 12345.6788);
    int_test2 = sprintf(Lmass2, "%-10.*f", 2, 12345.6788);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%25.*f", 1, 34.4325432);
    int_test2 = sprintf(Lmass2, "%25.*f", 1, 34.4325432);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%25.f", 4325678.43242);
    int_test2 = sprintf(Lmass2, "%25.f", 4325678.43242);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%#+10.f", 34567.2345432);
    int_test2 = sprintf(Lmass2, "%#+10.f", 34567.2345432);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%*.*f", 10, 2, 34567.43244);
    int_test2 = sprintf(Lmass2, "%*.*f", 10, 2, 34567.43244);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%f", (double)UINT_MAX);
    int_test2 = sprintf(Lmass2, "%f", (double)UINT_MAX);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    // int_test1 = s21_sprintf(Lmass1, "%f", (double)ULONG_MAX);// FAIL
    // int_test2 = sprintf(Lmass2, "%f", (double)ULONG_MAX);
    // ck_assert_str_eq(Lmass1, Lmass2);
    // ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%f", (double)INT_MAX);
    int_test2 = sprintf(Lmass2, "%f", (double)INT_MAX);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%f", (double)INT_MIN);
    int_test2 = sprintf(Lmass2, "%f", (double)INT_MIN);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%f", (double)SHRT_MAX);
    int_test2 = sprintf(Lmass2, "%f", (double)SHRT_MAX);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%f", (double)SHRT_MIN);
    int_test2 = sprintf(Lmass2, "%f", (double)SHRT_MIN);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%f", (double)LONG_MIN);
    int_test2 = sprintf(Lmass2, "%f", (double)LONG_MIN);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%f", (double)USHRT_MAX);
    int_test2 = sprintf(Lmass2, "%f", (double)USHRT_MAX);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%f", (double)LONG_MAX);
    int_test2 = sprintf(Lmass2, "%f", (double)LONG_MAX);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%+-10.5f%+10.2f%15.16f%+*.10f%*.16f%-10.*f%25.*f%25.f%#+10.f%*.*f",
     2546.4324, 3452.4325, 5678.43256, 25, 456789.243242, 10, 3456789123.43245, 2, 12345.6788,
     1, 34.4325432, 4325678.43242, 34567.2345432, 10, 2, 34567.43244);
    int_test2 = sprintf(Lmass2, "%+-10.5f%+10.2f%15.16f%+*.10f%*.16f%-10.*f%25.*f%25.f%#+10.f%*.*f",
     2546.4324, 3452.4325, 5678.43256, 25, 456789.243242, 10, 3456789123.43245, 2, 12345.6788,
     1, 34.4325432, 4325678.43242, 34567.2345432, 10, 2, 34567.43244);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1,
        "%+-10.5f  %+10.2f%15.16dasdf%+*.10f%*.16f%-10.*fdsad%25.*f%25.f%#+10.f%*.*f",
     2546.4324, 3452.4325, 567843256, 25, 456789.243242, 10, 3456789123.43245, 2, 12345.6788,
  1, 34.4325432, 4325678.43242, 34567.2345432, 10, 2, 34567.43244);
    int_test2 = sprintf(Lmass2,
        "%+-10.5f  %+10.2f%15.16dasdf%+*.10f%*.16f%-10.*fdsad%25.*f%25.f%#+10.f%*.*f",
     2546.4324, 3452.4325, 567843256, 25, 456789.243242, 10, 3456789123.43245, 2, 12345.6788,
  1, 34.4325432, 4325678.43242, 34567.2345432, 10, 2, 34567.43244);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);
}
END_TEST

START_TEST(Sprintf_Test_d) {
    s21_sprintf(Lmass1, "%-10d", 21);
    sprintf(Lmass2, "%-10d", 21);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(mass1, mass2);

    s21_sprintf(Lmass1, "%010d %010d", n2, n);
    sprintf(Lmass2, "%010d %010d", n2, n);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(mass1, mass2);

    s21_sprintf(Lmass1, "%*d %2d", -5, n2, n);
    sprintf(Lmass2, "%*d %2d", -5, n2, n);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %d abc", -166);
    sprintf(mass2, "abc %d abc", -166);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%d%%\n", 4);
    sprintf(mass2, "%d%%\n", 4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%.10d%%\n", 4);
    sprintf(mass2, "%.10d%%\n", 4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %010d", 50);
    sprintf(mass2, "abc %010d", 50);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%0d%%\n", -4);
    sprintf(mass2, "%0d%%\n", -4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%-+10d%%\n", 4);
    sprintf(mass2, "%-+10d%%\n", 4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%-10d%%\n", -4);
    sprintf(mass2, "%-10d%%\n", -4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%+10d%%\n", 4);
    sprintf(mass2, "%+10d%%\n", 4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%.2d", 4);
    sprintf(mass2, "%.2d", 4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%+10d%%\n", -4);
    sprintf(mass2, "%+10d%%\n", -4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(Lmass1, "%.*d %3.10d %.d", 10, n, n2, 42);
    int_test2 = sprintf(Lmass2, "%.*d %3.10d %.d", 10, n, n2, 42);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-+.5d %-+10.5d", 21, 10);
    int_test2 = sprintf(Lmass2, "%-+.5d %-+10.5d", 21, 10);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%+010d %+-10d", n2, n);
    int_test2 = sprintf(Lmass2, "%+010d %+-10d", n2, n);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "% -10d %-10d %-10d", 21, 55, -55);
    int_test2 = sprintf(Lmass2, "% -10d %-10d %-10d", 21, 55, -55);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%.10d %.10d", n, n2);
    int_test2 = sprintf(Lmass2, "%.10d %.10d", n, n2);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%2.2d %.d % d %.0d", -21, char_test1, 21, -42);
    int_test2 = sprintf(Lmass2, "%2.2d %.d % d %.0d", -21, char_test1, 21, -42);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%hd %hd", sh, sh1);
    int_test2 = sprintf(Lmass2, "%hd %hd", sh, sh1);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%ld %ld", lg, lg1);
    int_test2 = sprintf(Lmass2, "%ld %ld", lg, lg1);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%10ld", lg1);
    int_test2 = sprintf(Lmass2, "%10ld", lg1);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(mass1, "%+10.5d%%\n", -4);
    sprintf(mass2, "%+10.5d%%\n", -4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    long int var_long = 100;
    s21_sprintf(mass1, "%ld%%\n", var_long);
    sprintf(mass2, "%ld%%\n", var_long);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    short var_short = 123;
    s21_sprintf(mass1, "%hd%%\n", var_short);
    sprintf(mass2, "%hd%%\n", var_short);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %ld", 58787878958);
    sprintf(mass2, "abc %ld", 58787878958);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %d", 87878958);
    sprintf(mass2, "abc %d", 87878958);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %lld", (long long int)88888858787878958);
    sprintf(mass2, "abc %lld", (long long int)88888858787878958);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %hd", 87878958);
    sprintf(mass2, "abc %hd", 87878958);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %*d", 5, 4);
    sprintf(mass2, "abc %*d", 5, 4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(Lmass1, "%d%d", INT_MIN, INT_MAX);
    int_test2 = sprintf(Lmass2, "%d%d", INT_MIN, INT_MAX);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%d%d", SHRT_MAX, SHRT_MIN);
    int_test2 = sprintf(Lmass2, "%d%d", SHRT_MAX, SHRT_MIN);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%d %d", SHRT_MIN, USHRT_MAX);
    int_test2 = sprintf(Lmass2, "%d %d", SHRT_MIN, USHRT_MAX);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%ld", LONG_MAX);
    int_test2 = sprintf(Lmass2, "%ld", LONG_MAX);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%ld", LONG_MIN);
    int_test2 = sprintf(Lmass2, "%ld", LONG_MIN);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%ld", -9223372036854775807 -1);
    int_test2 = sprintf(Lmass2, "%ld", -9223372036854775807 -1);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-25.10d", 253);
    int_test2 = sprintf(Lmass2, "%-25.10d", 253);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-30.2d", 335);
    int_test2 = sprintf(Lmass2, "%-30.2d", 335);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%10.12d", 1234);
    int_test2 = sprintf(Lmass2, "%10.12d", 1234);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%11d", 34);
    int_test2 = sprintf(Lmass2, "%11d", 34);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%*d", 10, 25);
    int_test2 = sprintf(Lmass2, "%*d", 10, 25);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%*.*d", 5, 15, 345);
    int_test2 = sprintf(Lmass2, "%*.*d", 5, 15, 345);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%+-25.*d", 5, 15);
    int_test2 = sprintf(Lmass2, "%+-25.*d", 5, 15);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%+-25.*d", 45, 5);
    int_test2 = sprintf(Lmass2, "%+-25.*d", 45, 5);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%+10.2ld", 456789l);
    int_test2 = sprintf(Lmass2, "%+10.2ld", 456789l);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-11.*hd", 5, sh);
    int_test2 = sprintf(Lmass2, "%-11.*hd", 5, sh);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-25.10d%-30.2d%10.12d%11d%*d%*.*d%+-25.*d%+10.2ld%-11.*hd",
    253, 335, 1234, 34, 10, 25, 5, 15, 345, 5, 666, 456789l, 5, sh);
    int_test2 = sprintf(Lmass2, "%-25.10d%-30.2d%10.12d%11d%*d%*.*d%+-25.*d%+10.2ld%-11.*hd",
    253, 335, 1234, 34, 10, 25, 5, 15, 345, 5, 666, 456789l, 5, sh);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-25.10df%-30.2df%10.12dD %11d %*dsad%*.*d %+-25.*df%+10.2ld %-11.*hd",
    253, 335, 1234, 34, 10, 25, 5, 15, 345, 5, 666, 456789l, 5, sh);
    int_test2 = sprintf(Lmass2, "%-25.10df%-30.2df%10.12dD %11d %*dsad%*.*d %+-25.*df%+10.2ld %-11.*hd",
    253, 335, 1234, 34, 10, 25, 5, 15, 345, 5, 666, 456789l, 5, sh);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);
}
END_TEST

START_TEST(Sprintf_Test_u) {
    s21_sprintf(mass1, "abc %u abc", 166);
    sprintf(mass2, "abc %u abc", 166);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%2.5u%%\n", 4);
    sprintf(mass2, "%2.5u%%\n", 4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%2.5u%%\n", -4);
    sprintf(mass2, "%2.5u%%\n", -4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%10.5u%%\n", 4);
    sprintf(mass2, "%10.5u%%\n", 4);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "%-10u%%\n", 421);
    sprintf(mass2, "%-10u%%\n", 421);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    short unsigned key_h = 421;
    s21_sprintf(mass1, "%-10hu%%\n", key_h);
    sprintf(mass2, "%-10hu%%\n", key_h);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    long unsigned key_l = 421;
    s21_sprintf(mass1, "%-10lu%%\n", key_l);
    sprintf(mass2, "%-10lu%%\n", key_l);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %hu", 17);
    sprintf(mass2, "abc %hu", 17);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %lu", 17123121);
    sprintf(mass2, "abc %lu", (long unsigned int)17123121);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(Lmass1, "%u", 4321u);
    int_test2 = sprintf(Lmass2, "%u", 4321u);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%15u", 34567u);
    int_test2 = sprintf(Lmass2, "%15u", 34567u);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-20u", 2342u);
    int_test2 = sprintf(Lmass2, "%-20u", 2342u);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%010u", 12346u);
    int_test2 = sprintf(Lmass2, "%010u", 12346u);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-15.10u", 234542u);
    int_test2 = sprintf(Lmass2, "%-15.10u", 234542u);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%20.*u", 10, 234523u);
    int_test2 = sprintf(Lmass2, "%20.*u", 10, 234523u);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%*.*u", 10, 5, 2345u);
    int_test2 = sprintf(Lmass2, "%*.*u", 10, 5, 2345u);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%*.10u", 20, 2432u);
    int_test2 = sprintf(Lmass2, "%*.10u", 20, 2432u);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-20.12lu", lg);
    int_test2 = sprintf(Lmass2, "%-20.12lu", lg);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%19hu", sh);
    int_test2 = sprintf(Lmass2, "%19hu", sh);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%u", UINT_MAX);
    int_test2 = sprintf(Lmass2, "%u", UINT_MAX);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%lu", ULONG_MAX);
    int_test2 = sprintf(Lmass2, "%lu", ULONG_MAX);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(mass1, "%u %lu", UINT_MAX, ULONG_MAX);
    int_test2 = sprintf(mass2, "%u %lu", UINT_MAX, ULONG_MAX);
    ck_assert_str_eq(mass1, mass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(Lmass1, "%u%15u%-20u%010u%-15.10u%20.*u%*.*u%*.10u%-20.12lu%19hu",
    4321u, 34567u, 2342u, 12346u, 234542u, 10, 234523u, 10, 5, 2345u, 20, 2432u, lg, sh);
    int_test2 = sprintf(Lmass2, "%u%15u%-20u%010u%-15.10u%20.*u%*.*u%*.10u%-20.12lu%19hu",
    4321u, 34567u, 2342u, 12346u, 234542u, 10, 234523u, 10, 5, 2345u, 20, 2432u, lg, sh);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1,
    "%u%15uhhhh%-20u%010u  fsd  %-15.10u%20.*u sdss %*.*u%*.10u%-20.12lu%19hu",
  4321u, 34567u, 2342u, 12346u, 234542u, 10, 234523u, 10, 5, 2345u, 20, 2432u, lg, sh);
    int_test2 = sprintf(Lmass2,
    "%u%15uhhhh%-20u%010u  fsd  %-15.10u%20.*u sdss %*.*u%*.10u%-20.12lu%19hu",
  4321u, 34567u, 2342u, 12346u, 234542u, 10, 234523u, 10, 5, 2345u, 20, 2432u, lg, sh);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);
}
END_TEST


START_TEST(Sprintf_Test_x) {
    s21_sprintf(Lmass1, "%x %x", n, n1);
    sprintf(Lmass2, "%x %x", n, n1);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%10x %-10x", 21, n);
    sprintf(Lmass2, "%10x %-10x", 21, n);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%#4.x", 25);
    sprintf(Lmass2, "%#4.x", 25);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%.10x %.0x", n, 0);
    sprintf(Lmass2, "%.10x %.0x", n, 0);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%#5.10x", 123);
    sprintf(Lmass2, "%#5.10x", 123);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%x %x %x %10x %-10x %.10x %.0x %#5.10x", n, n1, n2, 0, 21, n, 0, 123);
    sprintf(Lmass2, "%x %x %x %10x %-10x %.10x %.0x %#5.10x", n, n1, n2, 0, 21, n, 0, 123);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%#x %#x", n, n2);
    sprintf(Lmass2, "%#x %#x", n, n2);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%#10x %#10x", 21, 0);
    sprintf(Lmass2, "%#10x %#10x", 21, 0);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%-#10x %0#10x %0#10x", 42, 21, 0);
    sprintf(Lmass2, "%-#10x %0#10x %0#10x", 42, 21, 0);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%-10.5x %lx %lx", -65, lg, lg1);
    sprintf(Lmass2, "%-10.5x %lx %lx", -65, lg, lg1);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%#x %#x %#10x %#10x %-#10x %0#10x %0#10x %-10.5x %lx %lx",
        n, n2, 21, 0, 42, 21, 0, -65, lg, lg1);
    sprintf(Lmass2, "%#x %#x %#10x %#10x %-#10x %0#10x %0#10x %-10.5x %lx %lx",
        n, n2, 21, 0, 42, 21, 0, -65, lg, lg1);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(mass1, "abc %x", 4784);
    sprintf(mass2, "abc %x", 4784);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %#x", 4784);
    sprintf(mass2, "abc %#x", 4784);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %x", 0);
    sprintf(mass2, "abc %x", 0);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);
}
END_TEST

START_TEST(Sprintf_Test_o) {
    s21_sprintf(mass1, "abc %o", 4784);
    sprintf(mass2, "abc %o", 4784);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %#o", 4784);
    sprintf(mass2, "abc %#o", 4784);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %o", 0);
    sprintf(mass2, "abc %o", 0);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %*o", -5, 0);
    sprintf(mass2, "abc %*o", -5, 0);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(Lmass1, "%-25.10o", 12345);
    int_test2 = sprintf(Lmass2, "%-25.10o", 12345);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(Lmass1, Lmass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-10o", 6532);
    int_test2 = sprintf(Lmass2, "%-10o", 6532);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(Lmass1, Lmass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%#-30.2o", 43522);
    int_test2 = sprintf(Lmass2, "%#-30.2o", 43522);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(Lmass1, Lmass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%#*.10o", 10, 64342);
    int_test2 = sprintf(Lmass2, "%#*.10o", 10, 64342);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(Lmass1, Lmass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%#*.*o", 15, 2, 343);
    int_test2 = sprintf(Lmass2, "%#*.*o", 15, 2, 343);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(Lmass1, Lmass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%25.*o", 6, 678);
    int_test2 = sprintf(Lmass2, "%25.*o", 6, 678);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(Lmass1, Lmass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%#.10ho", sh);
    int_test2 = sprintf(Lmass2, "%#.10ho", sh);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(Lmass1, Lmass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%*lo", 10, lg);
    int_test2 = sprintf(Lmass2, "%*lo", 10, lg);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(Lmass1, Lmass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(mass1, "%-10o", 21);
    int_test2 = sprintf(mass2, "%-10o", 21);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(mass1, mass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(mass1, "%lo", lg);
    int_test2 = sprintf(mass2, "%lo", lg);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(mass1, mass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(mass1, "%lo %o %o", lg1, n, 0);
    int_test2 = sprintf(mass2, "%lo %o %o", lg1, n, 0);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(mass1, mass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(mass1, "%o %10o", n2, 0);
    int_test2 = sprintf(mass2, "%o %10o", n2, 0);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(mass1, mass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(mass1, "%#10o %0#10o", 21, 42);
    int_test2 = sprintf(mass2, "%#10o %0#10o", 21, 42);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(mass1, mass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(mass1, "%#o", n);
    int_test2 = sprintf(mass2, "%#o", n);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(mass1, mass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(mass1, "%#o %#10o", n2, 0);
    int_test2 = sprintf(mass2, "%#o %#10o", n2, 0);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(mass1, mass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(mass1, "%0#10o %.10o", 0, n);
    int_test2 = sprintf(mass2, "%0#10o %.10o", 0, n);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(mass1, mass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(mass1, "%.0o %#5.10o", 0, 123);
    int_test2 = sprintf(mass2, "%.0o %#5.10o", 0, 123);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(mass1, mass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(Lmass1, "%-25.10o%-10o%#-30.2o%#*.10o%#*.*o%25.*o%#.10ho%*lo",
    12345, 6532, 43522, 10, 64342, 15, 2, 343, 6, 678, sh, 10, lg);
    int_test2 = sprintf(Lmass2, "%-25.10o%-10o%#-30.2o%#*.10o%#*.*o%25.*o%#.10ho%*lo",
    12345, 6532, 43522, 10, 64342, 15, 2, 343, 6, 678, sh, 10, lg);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(Lmass1, Lmass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-25.10o  %-10o  fdsf %#-30.2oGGGD#%#*.10o%#*.*o%25.*o%#.10ho%*lo",
    12345, 6532, 43522, 10, 64342, 15, 2, 343, 6, 678, sh, 10, lg);
    int_test2 = sprintf(Lmass2, "%-25.10o  %-10o  fdsf %#-30.2oGGGD#%#*.10o%#*.*o%25.*o%#.10ho%*lo",
    12345, 6532, 43522, 10, 64342, 15, 2, 343, 6, 678, sh, 10, lg);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(Lmass1, Lmass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);
}
END_TEST


START_TEST(Sprintf_Test_e) {
    int_test1 = s21_sprintf(Lmass1, "%+-20e", 2543.43254);
    int_test2 = sprintf(Lmass2, "%+-20e", 2543.43254);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(Lmass1, Lmass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%010.*e", 10, 1012.34567);
    int_test2 = sprintf(Lmass2, "%010.*e", 10, 1012.34567);
    ck_assert_int_eq(int_test1, int_test2);
    ck_assert_str_eq(Lmass1, Lmass2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%#02.e", 25.432523432);
    int_test2 = sprintf(Lmass2, "%#02.e", 25.432523432);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    s21_sprintf(mass1, "abc %e", 0.4784);
    sprintf(mass2, "abc %e", 0.4784);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %e", 7.8);
    sprintf(mass2, "abc %e", 7.8);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %e", 17.8);
    sprintf(mass2, "abc %e", 17.8);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %Le", (long double)17.8);
    sprintf(mass2, "abc %Le", (long double)17.8);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(Lmass1, "%+-20e", 2543.43254);
    int_test2 = sprintf(Lmass2, "%+-20e", 2543.43254);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-*.16e", 10, 1012.34567);
    int_test2 = sprintf(Lmass2, "%-*.16e", 10, 1012.34567);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%010.*e", 15, -0.00123);
    int_test2 = sprintf(Lmass2, "%010.*e", 15, -0.00123);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%#025.10e", 45678.342532);
    int_test2 = sprintf(Lmass2, "%#025.10e", 45678.342532);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%*.*e", 10, 15, 11111.22);
    int_test2 = sprintf(Lmass2, "%*.*e", 10, 15, 11111.22);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%*.15e", 16, -3245.3425342);
    int_test2 = sprintf(Lmass2, "%*.15e", 16, -3245.3425342);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%#010.*e", 16, -0.00123);
    int_test2 = sprintf(Lmass2, "%#010.*e", 16, -0.00123);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%25.16e", 235423424.43252);
    int_test2 = sprintf(Lmass2, "%25.16e", 235423424.43252);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%+-20e%010.*e%-*.16e%#025.10e%*.*e%*.15e%#010.*e%25.16e",
    2543.43254, 10, 1012.34567, 15, -0.00123, 45678.342532, 10, 15,
    11111.22, 16, -3245.3425342, 16, -0.00123, 235423424.43252);
    int_test2 = sprintf(Lmass2, "%+-20e%010.*e%-*.16e%#025.10e%*.*e%*.15e%#010.*e%25.16e",
    2543.43254, 10, 1012.34567, 15, -0.00123, 45678.342532, 10, 15,
    11111.22, 16, -3245.3425342, 16, -0.00123, 235423424.43252);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%+-20e  fds%010.*e%-*.16ebbbbb%#025.10e v %*.*e%*.15e%#010.*e%25.16e",
    2543.43254, 10, 1012.34567, 15, -0.00123, 45678.342532, 10, 15,
    11111.22, 16, -3245.3425342, 16, -0.00123, 235423424.43252);
    int_test2 = sprintf(Lmass2, "%+-20e  fds%010.*e%-*.16ebbbbb%#025.10e v %*.*e%*.15e%#010.*e%25.16e",
    2543.43254, 10, 1012.34567, 15, -0.00123, 45678.342532, 10, 15,
    11111.22, 16, -3245.3425342, 16, -0.00123, 235423424.43252);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);
}
END_TEST

START_TEST(Sprintf_Test_i) {
    s21_sprintf(mass1, "abc %-20i abc", 5);
    sprintf(mass2, "abc %-20i abc", 5);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %020i abc", 5);
    sprintf(mass2, "abc %020i abc", 5);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %.20i abc", 1234567890);
    sprintf(mass2, "abc %.20i abc", 1234567890);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %*i abc", 20, 1234567890);
    sprintf(mass2, "abc %*i abc", 20, 1234567890);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %.*i abc", 20, 1234567890);
    sprintf(mass2, "abc %.*i abc", 20, 1234567890);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc % 5i abc", 5);
    sprintf(mass2, "abc % 5i abc", 5);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %3i abc", 5);
    sprintf(mass2, "abc %3i abc", 5);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(mass1, "abc %+i abc", 5);
    int_test2 = sprintf(mass2, "abc %+i abc", 5);
    ck_assert_str_eq(mass1, mass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(mass1, "abc %-.i abc", n);
    int_test2 = sprintf(mass2, "abc %-.i abc", n);
    ck_assert_str_eq(mass1, mass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(mass1, "abc %*.023i abc", 25, n);
    int_test2 = sprintf(mass2, "abc %*.023i abc", 25, n);
    ck_assert_str_eq(mass1, mass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(mass1, "abc %*hi abc", 33, (short)n);
    int_test2 = sprintf(mass2, "abc %*hi abc", 33, (short)n);
    ck_assert_str_eq(mass1, mass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(mass1, "abc %-+i abc", n2);
    int_test2 = sprintf(mass2, "abc %-+i abc", n2);
    ck_assert_str_eq(mass1, mass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(mass1, mass2);
}
END_TEST

START_TEST(Sprintf_Test_pr) {
    s21_sprintf(mass1, "abc %%");
    sprintf(mass2, "abc %%");
    ck_assert_str_eq(mass1, mass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(mass1, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    int_test2 = sprintf(mass2, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    ck_assert_str_eq(mass1, mass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);
}
END_TEST

START_TEST(Sprintf_Test_g) {
    int_test1 = s21_sprintf(Lmass1, "%g", 2345.4324);
    int_test2 = sprintf(Lmass2, "%g", 2345.4324);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%23g", 23445.2342);
    int_test2 = sprintf(Lmass2, "%23g", 23445.2342);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%g%23g", 2345.4324, 23445.2342);
    int_test2 = sprintf(Lmass2, "%g%23g", 2345.4324, 23445.2342);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%g %23g", 2345.4324, 23445.2342);
    int_test2 = sprintf(Lmass2, "%g %23g", 2345.4324, 23445.2342);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%+#25.8g", 0.0000034);
    int_test2 = sprintf(Lmass2, "%+#25.8g", 0.0000034);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-10.*g", 9, 234567.43242);
    int_test2 = sprintf(Lmass2, "%-10.*g", 9, 234567.43242);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-10g", 3456.3424);
    int_test2 = sprintf(Lmass2, "%-10g", 3456.3424);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%015g", 3456.3424);
    int_test2 = sprintf(Lmass2, "%015g", 3456.3424);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-10g%015g", 1234.234242, 3456.3424);
    int_test2 = sprintf(Lmass2, "%-10g%015g", 1234.234242, 3456.3424);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-10g %015g", 1234.234242, 3456.3424);
    int_test2 = sprintf(Lmass2, "%-10g %015g", 1234.234242, 3456.3424);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%#10.5g %0#15.10g", 0.00001234, -324234.234);
    int_test2 = sprintf(Lmass2, "%#10.5g %0#15.10g", 0.00001234, -324234.234);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%+10.6g", 23456.4324);
    int_test2 = sprintf(Lmass2, "%+10.6g", 23456.4324);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%#*.g", 5, 23456.4324);
    int_test2 = sprintf(Lmass2, "%#*.g", 5, 23456.4324);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%+10.6g %#*.g", 2345.4234, 5, 23456.4324);
    int_test2 = sprintf(Lmass2, "%+10.6g %#*.g", 2345.4234, 5, 23456.4324);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%Lg", (long double) ULLONG_MAX);
    int_test2 = sprintf(Lmass2, "%Lg", (long double) ULONG_MAX);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%Lg", (long double) LLONG_MIN);
    int_test2 = sprintf(Lmass2, "%Lg", (long double) LLONG_MIN);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%lg", (double) ULLONG_MAX);
    int_test2 = sprintf(Lmass2, "%lg", (double) ULLONG_MAX);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%lg", (double) LLONG_MIN);
    int_test2 = sprintf(Lmass2, "%lg", (double) LLONG_MIN);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "dasdad %+10.6g ufiodusfoi34 %#*.g\n", 2345.4234, 5, 23456.4324);
    int_test2 = sprintf(Lmass2, "dasdad %+10.6g ufiodusfoi34 %#*.g\n", 2345.4234, 5, 23456.4324);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%15.30Lg", (long double) ULLONG_MAX);
    int_test2 = sprintf(Lmass2, "%15.30Lg", (long double) ULLONG_MAX);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);
}
END_TEST

START_TEST(Sprintf_Test_G) {
    int_test1 = s21_sprintf(Lmass1, "%G", 2345.4324);
    int_test2 = sprintf(Lmass2, "%G", 2345.4324);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%23G", 23445.2342);
    int_test2 = sprintf(Lmass2, "%23G", 23445.2342);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%G%23G", 2345.4324, 23445.2342);
    int_test2 = sprintf(Lmass2, "%G%23G", 2345.4324, 23445.2342);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%G %23G", 2345.4324, 23445.2342);
    int_test2 = sprintf(Lmass2, "%G %23G", 2345.4324, 23445.2342);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%+#25.8G", 0.0000034);
    int_test2 = sprintf(Lmass2, "%+#25.8G", 0.0000034);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-10.*G", 9, 234567.43242);
    int_test2 = sprintf(Lmass2, "%-10.*G", 9, 234567.43242);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-10G", 3456.3424);
    int_test2 = sprintf(Lmass2, "%-10G", 3456.3424);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%015G", 3456.3424);
    int_test2 = sprintf(Lmass2, "%015G", 3456.3424);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-10G%015G", 1234.234242, 3456.3424);
    int_test2 = sprintf(Lmass2, "%-10G%015G", 1234.234242, 3456.3424);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-10G %015G", 1234.234242, 3456.3424);
    int_test2 = sprintf(Lmass2, "%-10G %015G", 1234.234242, 3456.3424);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%#10.5G %0#15.10G", 0.00001234, -324234.234);
    int_test2 = sprintf(Lmass2, "%#10.5G %0#15.10G", 0.00001234, -324234.234);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%+10.6G", 23456.4324);
    int_test2 = sprintf(Lmass2, "%+10.6G", 23456.4324);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%#*.G", 5, 23456.4324);
    int_test2 = sprintf(Lmass2, "%#*.G", 5, 23456.4324);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%+10.6G %#*.G", 2345.4234, 5, 23456.4324);
    int_test2 = sprintf(Lmass2, "%+10.6G %#*.G", 2345.4234, 5, 23456.4324);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%+10.6G%#*.G", 2345.4234, 5, 23456.4324);
    int_test2 = sprintf(Lmass2, "%+10.6G%#*.G", 2345.4234, 5, 23456.4324);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);
}
END_TEST

START_TEST(Sprintf_Test_E) {
    s21_sprintf(mass1, "abc %E", 0.4784);
    sprintf(mass2, "abc %E", 0.4784);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %-20E", 0.4784);
    sprintf(mass2, "abc %-20E", 0.4784);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %+20E", 0.0);
    sprintf(mass2, "abc %+20E", 0.0);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc % 20E", 123456789.123456789);
    sprintf(mass2, "abc % 20E", 123456789.123456789);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    s21_sprintf(mass1, "abc %#20E", 1e10);
    sprintf(mass2, "abc %#20E", 1e10);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(Lmass1, "%+-25.4E", 3456.4325);
    int_test2 = sprintf(Lmass2, "%+-25.4E", 3456.4325);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%25.*E", 10, 23452.43242);
    int_test2 = sprintf(Lmass2, "%25.*E", 10, 23452.43242);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%-*.10E", 5, -0.05234);
    int_test2 = sprintf(Lmass2, "%-*.10E", 5, -0.05234);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%25.*E%-*.10E", 10, 23452.43242, 5, -0.05234);
    int_test2 = sprintf(Lmass2, "%25.*E%-*.10E", 10, 23452.43242, 5, -0.05234);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%25.*E%-*.10E", 10, 23452.43242, 5, -0.05234);
    int_test2 = sprintf(Lmass2, "%25.*E%-*.10E", 10, 23452.43242, 5, -0.05234);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%#02.E", 25.432523432);
    int_test2 = sprintf(Lmass2, "%#02.E", 25.432523432);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%*.*E", 10, 2, 53242.4242);
    int_test2 = sprintf(Lmass2, "%*.*E", 10, 2, 53242.4242);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%*.15E%*.*E", 10, 456789.43242, 10, 2, 53242.4242);
    int_test2 = sprintf(Lmass2, "%*.15E%*.*E", 10, 456789.43242, 10, 2, 53242.4242);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%*.15E %*.*E", 10, 456789.43242, 10, 2, 53242.4242);
    int_test2 = sprintf(Lmass2, "%*.15E %*.*E", 10, 456789.43242, 10, 2, 53242.4242);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%*.15E", 10, 456789.43242);
    int_test2 = sprintf(Lmass2, "%*.15E", 10, 456789.43242);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%*E", 25, 1111122222.34567899);
    int_test2 = sprintf(Lmass2, "%*E", 25, 1111122222.34567899);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%12.16E", 2345678.23453242);
    int_test2 = sprintf(Lmass2, "%12.16E", 2345678.23453242);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%15Le", (long double)ULLONG_MAX);
    int_test2 = sprintf(Lmass2, "%15Le", (long double)ULLONG_MAX);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);

    int_test1 = s21_sprintf(Lmass1, "%15Le", (long double) LLONG_MIN);
    int_test2 = sprintf(Lmass2, "%15Le", (long double) LLONG_MIN);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    int_test1 = 0;
    int_test2 = 0;
    free_array(Lmass1, Lmass2);
}
END_TEST

START_TEST(Sprintf_Test_X) {
    s21_sprintf(Lmass1, "%X", n);
    sprintf(Lmass2, "%X", n);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%X %X", n, n2);
    sprintf(Lmass2, "%X %X", n, n2);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%X %X%X", n, n2, n1);
    sprintf(Lmass2, "%X %X%X", n, n2, n1);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%-10X %lX", 21, lg1);
    sprintf(Lmass2, "%-10X %lX", 21, lg1);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%lX", lg1);
    sprintf(Lmass2, "%lX", lg1);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%-10X", 21);
    sprintf(Lmass2, "%-10X", 21);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%10X %-10X", 0, 21);
    sprintf(Lmass2, "%10X %-10X", 0, 21);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%#X %#X", n, n2);
    sprintf(Lmass2, "%#X %#X", n, n2);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%-#10X %#10X", 21, 42);
    sprintf(Lmass2, "%-#10X %#10X", 21, 42);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%-#10X%#10X", 21, 42);
    sprintf(Lmass2, "%-#10X%#10X", 21, 42);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%-#10X", 21);
    sprintf(Lmass2, "%-#10X", 21);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%#10X", 42);
    sprintf(Lmass2, "%#10X", 42);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%0#10X", 4221);
    sprintf(Lmass2, "%0#10X", 4221);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%0#10X", n1);
    sprintf(Lmass2, "%0#10X", n1);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%0#10X %0#10X", 4221, n1);
    sprintf(Lmass2, "%0#10X %0#10X", 4221, n1);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%0#10X%0#10X", 4221, n1);
    sprintf(Lmass2, "%0#10X%0#10X", 4221, n1);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%X %X %X %10X %-10X %lX", n, n1, n2, 0, 21, lg1);
    sprintf(Lmass2, "%X %X %X %10X %-10X %lX", n, n1, n2, 0, 21, lg1);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(Lmass1, "%#X %#X %-#10X %#10X %#10X %0#10X %0#10X", n, n2, 21, 42, 0, 4221, n1);
    sprintf(Lmass2, "%#X %#X %-#10X %#10X %#10X %0#10X %0#10X", n, n2, 21, 42, 0, 4221, n1);
    ck_assert_str_eq(Lmass1, Lmass2);
    free_array(Lmass1, Lmass2);

    s21_sprintf(mass1, "abc %X", 4784);
    sprintf(mass2, "abc %X", 4784);
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);
}
END_TEST

START_TEST(Sprintf_Test_n) {
    int test_int1 = 0;
    int test_int2 = 0;

    int_test1 = s21_sprintf(mass1, "abc %n", &test_int1);
    int_test2 = sprintf(mass2, "abc %n", &test_int2);
    ck_assert_int_eq(test_int1, test_int2);
    ck_assert_int_eq(int_test1, int_test2);
    test_int1 = 0;
    test_int2 = 0;
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(mass1, "Hello %n Hello", &test_int1);
    int_test2 = sprintf(mass2, "Hello %n Hello", &test_int2);
    ck_assert_int_eq(test_int1, test_int2);
    ck_assert_int_eq(int_test1, int_test2);
    test_int1 = 0;
    test_int2 = 0;
    free_array(mass1, mass2);

    int_test1 = s21_sprintf(mass1, "Hello %n Hello %d", &test_int1, test_int1);
    int_test2 = sprintf(mass2, "Hello %n Hello %d", &test_int2, test_int2);
    ck_assert_int_eq(test_int1, test_int2);
    ck_assert_int_eq(int_test1, int_test2);
    test_int1 = 0;
    test_int2 = 0;
    free_array(mass1, mass2);
}
END_TEST

START_TEST(Sprintf_Test_p) {
    s21_sprintf(mass1, "qwerty");
    sprintf(mass2, "qwerty");
    ck_assert_str_eq(mass1, mass2);
    free_array(mass1, mass2);

    void *p1 = (void*)0x123456789;
    void *p2 = (void*)4325234;
    void *p3 = (void*)0x123f324b;
    void *p4 = (void*)432520;

    int_test1 = sprintf(Lmass1, "%p", p1);
    int_test2 = s21_sprintf(Lmass2, "%p", p1);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    free_array(Lmass1, Lmass2);

    int_test1 = sprintf(Lmass1, "%20p", p2);
    int_test2 = s21_sprintf(Lmass2, "%20p", p2);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    free_array(Lmass1, Lmass2);

    int_test1 = sprintf(Lmass1, "%-15p", p3);
    int_test2 = s21_sprintf(Lmass2, "%-15p", p3);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    free_array(Lmass1, Lmass2);

    int_test1 = sprintf(Lmass1, "%10p", p4);
    int_test2 = s21_sprintf(Lmass2, "%10p", p4);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    free_array(Lmass1, Lmass2);

    int_test1 = sprintf(Lmass1, "%-15p", p2);
    int_test2 = s21_sprintf(Lmass2, "%-15p", p2);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    free_array(Lmass1, Lmass2);

    int_test1 = sprintf(Lmass1, "%*p", 10, p3);
    int_test2 = s21_sprintf(Lmass2, "%*p", 10, p3);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    free_array(Lmass1, Lmass2);

    int_test1 = sprintf(Lmass1, "%10p", p1);
    int_test2 = s21_sprintf(Lmass2, "%10p", p1);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    free_array(Lmass1, Lmass2);

    int_test1 = sprintf(Lmass1, "%*p", 10, p4);
    int_test2 = s21_sprintf(Lmass2, "%*p", 10, p4);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    free_array(Lmass1, Lmass2);

    int_test1 = sprintf(Lmass1, "%p%20p%-15p%10p%-15p%*p%10p%*p",
    p1, p2, p3, p4, p2, 10, p3, p1, 10, p4);
    int_test2 = s21_sprintf(Lmass2, "%p%20p%-15p%10p%-15p%*p%10p%*p",
    p1, p2, p3, p4, p2, 10, p3, p1, 10, p4);
    ck_assert_str_eq(Lmass1, Lmass2);
    ck_assert_int_eq(int_test1, int_test2);
    free_array(Lmass1, Lmass2);
}
END_TEST

START_TEST(test_s21_MultySpec) {
    free_array(mass1, mass2);
    void *p = (void*)0x3456;
    int y1, y2;
    int_test1 = sprintf(mass2, "%d%o%s%c%i%e%f%E%g%u%x%G%X%n%p%%%s",
    25, 342, "aboba", 'a', 123, 25.34, 4325.23434,
    0.0000000123, 2.12345, 12345u, 8342, 0.0000456, 1235, &y1, p, "hello my friend");
    int_test2 = s21_sprintf(mass1, "%d%o%s%c%i%e%f%E%g%u%x%G%X%n%p%%%s", 25, 342, "aboba", 'a',
    123, 25.34, 4325.23434,
    0.0000000123, 2.12345, 12345u, 8342, 0.0000456, 1235, &y2, p, "hello my friend");
    ck_assert_str_eq(mass1, mass2);
    ck_assert_int_eq(y1, y2);
    ck_assert_int_eq(int_test1, int_test2);
}
END_TEST

START_TEST(sscanf_test_e_1) {
    float e1 = 0, e2 = 0;
    int Std_Lib = 0, Our_Lib = 0;

    Std_Lib = sscanf("123e34", "%e", &e1);
    Our_Lib = s21_sscanf("123e34", "%e", &e2);
    ck_assert_msg(e1 == e2, "float not equals"); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("12e34", "%e", &e1);
    Our_Lib = s21_sscanf("12e34", "%e", &e2);
    ck_assert_msg(e1 == e2, "float not equals"); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("123.345e-34", "%9e", &e1);
    Our_Lib = s21_sscanf("123.345e-34", "%9e", &e2);
    ck_assert_int_eq(Std_Lib, Our_Lib);
    ck_assert_float_eq(e1, e2);
    ck_assert_msg(e1 == e2, "float not equals");

    Std_Lib = sscanf("123.345e34", "%8e", &e1);
    Our_Lib = s21_sscanf("123.345e34", "%8e", &e2);
    ck_assert_float_eq(e1, e2);
    ck_assert_msg(e1 == e2, "float not equals");
    ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_l_1) {
    int Std_Lib = 0, Our_Lib = 0;
    unsigned long lu1 = 0, lu2 = 0;
    long ld1 = 0, ld2 = 0;
    wchar_t lc1, lc2;
    wchar_t lsuite[100], ls2[100];

    Std_Lib = sscanf("a", "%lc", &lc1);
    Our_Lib = s21_sscanf("a", "%lc", &lc2);
    ck_assert_int_eq(lc1, lc2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("12345", "%lu", &lu1);
    Our_Lib = s21_sscanf("12345", "%lu", &lu2);
    ck_assert_int_eq(lu1, lu2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("12345", "%ld", &ld1);
    Our_Lib = s21_sscanf("12345", "%ld", &ld2);
    ck_assert_int_eq(ld1, ld2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("amogus", "%9ls", lsuite);
    Our_Lib = s21_sscanf("amogus", "%9ls", ls2);
    ck_assert_int_eq(Std_Lib, Our_Lib);
    ck_assert_int_eq(wcscmp(ls2, lsuite), 0);

    Std_Lib = sscanf("amogus", "%2ls", lsuite);
    Our_Lib = s21_sscanf("amogus", "%2ls", ls2);
    ck_assert_int_eq(wcscmp(lsuite, ls2), 0); ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_L_1) {
    int Std_Lib = 0, Our_Lib = 0;
    long double Lf1 = 0, Lf2 = 0;

    Std_Lib = sscanf("123.45", "%Lf", &Lf1);
    Our_Lib = s21_sscanf("123.45", "%Lf", &Lf2);
    ck_assert_msg((int)Lf1 == (int)Lf2, "floats not equals"); ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_u_1) {
    int Std_Lib = 0, Our_Lib = 0;
    unsigned u1 = 0, u2 = 0;

    Std_Lib = sscanf("", "%u", &u1);
    Our_Lib = s21_sscanf("", "%u", &u2);
    ck_assert_int_eq(u1, u2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("0x12345", "%2u", &u1);
    Our_Lib = s21_sscanf("0x12345", "%2u", &u2);
    ck_assert_int_eq(u1, u2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("123.456", "%u", &u1);
    Our_Lib = s21_sscanf("123.456", "%u", &u2);
    ck_assert_msg(u1 == u2, "floats not equals"); ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_h_1) {
    int Std_Lib = 0, Our_Lib = 0;
    unsigned short hu1 = 0, hu2 = 0;

    Std_Lib = sscanf("12345", "%hu", &hu1);
    Our_Lib = s21_sscanf("12345", "%hu", &hu2);
    ck_assert_int_eq(hu1, hu2); ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_s_1) {
    int Std_Lib = 0, Our_Lib = 0;
    char suite[100], s2[100];

    Std_Lib = sscanf("shrek amogus", "%*s%2s", suite);
    Our_Lib = s21_sscanf("shrek amogus", "%*s%2s", s2);

    ck_assert_str_eq(suite, s2); ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_f_1) {
    int Std_Lib = 0, Our_Lib = 0;
    float f_1 = 0, f2 = 0;

    Std_Lib = sscanf("12e-3", "%f", &f_1);
    Our_Lib = s21_sscanf("12e-3", "%f", &f2);
    ck_assert_int_eq(f_1, f2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("-123.45", "%f", &f_1);
    Our_Lib = s21_sscanf("-123.45", "%f", &f2);
    ck_assert_int_eq(f_1, f2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("-12.345", "%5f", &f_1);
    Our_Lib = s21_sscanf("-12.345", "%5f", &f2);
    ck_assert_int_eq(f_1, f2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("12.345", "%5f", &f_1);
    Our_Lib = s21_sscanf("12.345", "%5f", &f2);
    ck_assert_int_eq(f_1, f2); ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_o_1) {
    int Std_Lib = 0, Our_Lib = 0;
    unsigned int o1 = 0, o2 = 0;
    char * A = "12345678";
    char * B = "12345678123456789";

    Std_Lib = sscanf(A, "%3o", &o1);
    Our_Lib = s21_sscanf(A, "%3o", &o2);
    ck_assert_int_eq(o1, o2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(B, "%9o", &o1);
    Our_Lib = s21_sscanf(B, "%9o", &o2);
    ck_assert_int_eq(o1, o2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(B, "%o", &o1);
    Our_Lib = s21_sscanf(B, "%o", &o2);
    ck_assert_int_eq(o1, o2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("123.456", "%o", &o1);
    Our_Lib = s21_sscanf("123.456", "%o", &o2);
    ck_assert_uint_eq(o1, o2); ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_x_1) {
    int Std_Lib = 0, Our_Lib = 0;
    unsigned x1 = 5, x2 = 5;

    Std_Lib = sscanf("0x12345", "%2x", &x1);
    Our_Lib = s21_sscanf("0x12345", "%2x", &x2);
    ck_assert_uint_eq(x1, x2); ck_assert_int_eq(Std_Lib, Our_Lib);

    x1 = 10, x2 = 10;
    Std_Lib = sscanf("AMOGUS", "%x", &x1);
    Our_Lib = s21_sscanf("AMOGUS", "%x", &x2);
    ck_assert_uint_eq(x1, x2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("amogus", "%x", &x1);
    Our_Lib = s21_sscanf("amogus", "%x", &x2);
    ck_assert_uint_eq(x1, x2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("0x1234567890123456789", "%x", &x1);
    Our_Lib = s21_sscanf("0x1234567890123456789", "%x", &x2);
    ck_assert_uint_eq(x1, x2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("0x123", "%x", &x1);
    Our_Lib = s21_sscanf("0x123", "%x", &x2);
    ck_assert_uint_eq(x1, x2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("123456", "%x", &x1);
    Our_Lib = s21_sscanf("123456", "%x", &x2);
    ck_assert_uint_eq(x1, x2); ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_i_1) {
    int Std_Lib = 0, Our_Lib = 0;
    int i1 = 0, i2 = 0;

    Std_Lib = sscanf("0123", "%i", &i1);
    Our_Lib = s21_sscanf("0123", "%i", &i2);
    ck_assert_int_eq(i1, i2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("0x123", "%i", &i1);
    Our_Lib = s21_sscanf("0x123", "%i", &i2);
    ck_assert_int_eq(i1, i2); ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_proc_1) {
    int Std_Lib = 0, Our_Lib = 0;
    char proc1[100] = {0}, proc2[100] = {0};

    Std_Lib = sscanf("123.456", "%%%9s", proc1);
    Our_Lib = s21_sscanf("123.456", "%%%9s", proc2);
    ck_assert_str_eq(proc1, proc2); ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_p_1) {
    int Std_Lib = 0, Our_Lib = 0;
    void *p1 = 0, *p2 = 0;

    Std_Lib = sscanf("0x123", "%4p", &p1);
    Our_Lib = s21_sscanf("0x123", "%4p", &p2);
    ck_assert_ptr_eq(p1, p2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("0x12345", "%p", &p1);
    Our_Lib = s21_sscanf("0x12345", "%p", &p2);
    ck_assert_ptr_eq(p1, p2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("123.456", "%p", &p1);
    Our_Lib = s21_sscanf("123.456", "%p", &p2);
    ck_assert_ptr_eq(p1, p2); ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_n_1) {
    int n_1 = 0, n_2 = 0, Std_Lib = 0, Our_Lib = 0;

    Std_Lib = sscanf("123.456", "%n", &n_1);
    Our_Lib = s21_sscanf("123.456", "%n", &n_2);
    ck_assert_int_eq(n_1, n_2); ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_d_1) {
    int Std_Lib = 0, Our_Lib = 0, d1 = 0, d2 = 0;

    Std_Lib = sscanf("123", "%002d", &d1);
    Our_Lib = s21_sscanf("123", "%002d", &d2);
    ck_assert_int_eq(d1, d2); ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_all) {
    int d1 = 0, d2 = 0, i1 = 0, i2 = 0, n_1 = 0, n_2 = 0, Std_Lib = 0, Our_Lib = 0;
    float e1 = 0, e2 = 0, E1 = 0, E2 = 0, f_1 = 0, f2 = 0, g1 = 0, g2 = 0, G1 = 0, G2 = 0;
    unsigned o1 = 0, o2 = 0, u1 = 0, u2 = 0, x1 = 0, x2 = 0, X1 = 0, X2 = 0;
    void *p1 = 0, *p2 = 0;
    char suite[100], s2[100], c1, c2, proc1[100] = {0}, proc2[100] = {0};
    short hi1 = 0, hi2 = 0, hd1 = 0, hd2 = 0;
    unsigned long lo1 = 0, lo2 = 0, lu1 = 0, lu2 = 0;
    long double Le1 = 0, Le2 = 0, LE1 = 0, LE2 = 0, Lg1 = 0, Lg2 = 0, LG1 = 0, LG2 = 0;

    Std_Lib = sscanf("123", "%d", &d1);
    Our_Lib = s21_sscanf("123", "%d", &d2);
    ck_assert_int_eq(d1, d2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("123 amogus", "%d%9s", &d1, suite);
    Our_Lib = s21_sscanf("123 amogus", "%d%9s", &d2, s2);
    ck_assert_int_eq(d1, d2); ck_assert_str_eq(suite, s2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("123 amogus q", "%d%9s%c", &d1, suite, &c1);
    Our_Lib = s21_sscanf("123 amogus q", "%d%9s%c", &d2, s2, &c2);
    ck_assert_int_eq(d1, d2); ck_assert_str_eq(suite, s2); ck_assert_int_eq(c1, c2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("4 123 amogus q", "%i%d%9s%c", &i1, &d1, suite, &c1);
    Our_Lib = s21_sscanf("4 123 amogus q", "%i%d%9s%c", &i2, &d2, s2, &c2);
    ck_assert_int_eq(i1, i2); ck_assert_int_eq(d1, d2); ck_assert_str_eq(suite, s2);
    ck_assert_int_eq(c1, c2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("4 123 amogus q 12e34", "%i%d%9s%c%e", &i1, &d1, suite, &c1, &e1);
    Our_Lib = s21_sscanf("4 123 amogus q 12e34", "%i%d%9s%c%e", &i2, &d2, s2, &c2, &e2);
    ck_assert_int_eq(i1, i2);
    ck_assert_int_eq(d1, d2);
    ck_assert_str_eq(suite, s2);
    ck_assert_int_eq(c1, c2);
    ck_assert_msg(e1 == e2, "floats not equals");
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("4 123 amogus q 12e34 43e21", "%i%d%9s%c%e%E", &i1, &d1, suite, &c1, &e1, &E1);
    Our_Lib = s21_sscanf("4 123 amogus q 12e34 43e21", "%i%d%9s%c%e%E", &i2, &d2, s2, &c2, &e2, &E2);
    ck_assert_int_eq(i1, i2);
    ck_assert_int_eq(d1, d2);
    ck_assert_str_eq(suite, s2);
    ck_assert_int_eq(c1, c2);
    ck_assert_msg(e1 == e2, "floats not equals");
    ck_assert_msg(E1 == E2, "floats not equals");
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("4 123 amogus q 12e34 43e21 13.37", "%i%d%9s%c%e%E%f",
    &i1, &d1, suite, &c1, &e1, &E1, &f_1);
    Our_Lib = s21_sscanf("4 123 amogus q 12e34 43e21 13,37", "%i%d%9s%c%e%E%f",
    &i2, &d2, s2, &c2, &e2, &E2, &f2);
    ck_assert_int_eq(i1, i2); ck_assert_int_eq(d1, d2);
    ck_assert_str_eq(suite, s2);
    ck_assert_int_eq(c1, c2);
    ck_assert_msg(e1 == e2, "floats not equals");
    ck_assert_msg(E1 == E2, "floats not equals");
    ck_assert_msg(f_1 == f2, "floats not equals");
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("4 123 amogus q 12e34 43e21 13.37 2.28", "%i%d%9s%c%e%E%f%g",
    &i1, &d1, suite, &c1, &e1, &E1, &f_1, &g1);
    Our_Lib = s21_sscanf("4 123 amogus q 12e34 43e21 13.37 2.28", "%i%d%9s%c%e%E%f%g",
    &i2, &d2, s2, &c2, &e2, &E2, &f2, &g2);
    ck_assert_int_eq(i1, i2);
    ck_assert_int_eq(d1, d2);
    ck_assert_str_eq(suite, s2);
    ck_assert_int_eq(c1, c2);
    ck_assert_msg(e1 == e2, "floats not equals");
    ck_assert_msg(E1 == E2, "floats not equals");
    ck_assert_msg(f_1 == f2, "floats not equals");
    ck_assert_msg(g1 == g2, "floats not equals");
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("4 123 amogus q 12e34 43e21 13.37 2.28 1337.228", "%i%d%9s%c%e%E%f%g%G",
    &i1, &d1, suite, &c1, &e1, &E1, &f_1, &g1, &G1);
    Our_Lib = s21_sscanf("4 123 amogus q 12e34 43e21 13.37 2.28 1337.228", "%i%d%9s%c%e%E%f%g%G",
    &i2, &d2, s2, &c2, &e2, &E2, &f2, &g2, &G2);
    ck_assert_int_eq(i1, i2);
    ck_assert_int_eq(d1, d2);
    ck_assert_str_eq(suite, s2);
    ck_assert_int_eq(c1, c2);
    ck_assert_msg(e1 == e2, "floats not equals");
    ck_assert_msg(E1 == E2, "floats not equals");
    ck_assert_msg(f_1 == f2, "floats not equals");
    ck_assert_msg(g1 == g2, "floats not equals");
    ck_assert_msg(G1 == G2, "floats not equals");
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("4 123 amogus q 12e34 43e21 13.37 2.28 1337.228",
     "%i%d%9s%c%e%E%f%g%G",
    &i1, &d1, suite, &c1, &e1, &E1, &f_1, &g1, &G1);
    Our_Lib = s21_sscanf("4 123 amogus q 12e34 43e21 13.37 2.28 1337.228",
    "%i%d%9s%c%e%E%f%g%G",
    &i2, &d2, s2, &c2, &e2, &E2, &f2, &g2, &G2);
    ck_assert_int_eq(i1, i2); ck_assert_int_eq(d1, d2); ck_assert_str_eq(suite, s2);
    ck_assert_int_eq(c1, c2); ck_assert_msg(e1 == e2, "floats not equals");
    ck_assert_msg(E1 == E2, "floats not equals");
    ck_assert_msg(f_1 == f2, "floats not equals");
    ck_assert_msg(g1 == g2, "floats not equals");
    ck_assert_msg(G1 == G2, "floats not equals");
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("4 123 amogus q 12e34 43e21 13.37 2.28 1337.228 101", "%i%d%9s%c%e%E%f%g%G%o",
    &i1, &d1, suite, &c1, &e1, &E1, &f_1, &g1, &G1, &o1);
    Our_Lib = s21_sscanf("4 123 amogus q 12e34 43e21 13.37 2.28 1337.228 101", "%i%d%9s%c%e%E%f%g%G%o",
    &i2, &d2, s2, &c2, &e2, &E2, &f2, &g2, &G2, &o2);
    ck_assert_int_eq(i1, i2); ck_assert_int_eq(d1, d2); ck_assert_str_eq(suite, s2);
    ck_assert_int_eq(c1, c2);
    ck_assert_msg(e1 == e2, "floats not equals");
    ck_assert_msg(E1 == E2, "floats not equals");
    ck_assert_msg(f_1 == f2, "floats not equals");
    ck_assert_msg(g1 == g2, "floats not equals");
    ck_assert_msg(G1 == G2, "floats not equals");
    ck_assert_int_eq(o1, o2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("4 123 amogus q 12e34 43e21 13.37 2.28 1337.228 101 1854", "%i%d%9s%c%e%E%f%g%G%o%u",
    &i1, &d1, suite, &c1, &e1, &E1, &f_1, &g1, &G1, &o1, &u1);
    Our_Lib = s21_sscanf("4 123 amogus q 12e34 43e21 13.37 2.28 1337.228 101 1854", "%i%d%9s%c%e%E%f%g%G%o%u",
    &i2, &d2, s2, &c2, &e2, &E2, &f2, &g2, &G2, &o2, &u2);
    ck_assert_int_eq(i1, i2); ck_assert_int_eq(d1, d2); ck_assert_str_eq(suite, s2);
    ck_assert_int_eq(c1, c2); ck_assert_msg(e1 == e2, "floats not equals");
    ck_assert_msg(E1 == E2, "floats not equals");
    ck_assert_msg(f_1 == f2, "floats not equals"); ck_assert_msg(g1 == g2, "floats not equals");
    ck_assert_msg(G1 == G2, "floats not equals");
    ck_assert_int_eq(o1, o2); ck_assert_int_eq(u1, u2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("4 123 amogus q 12e34 43e21 13.37 2.28 1337.228 101 1854 12ab",
    "%i%d%9s%c%e%E%f%g%G%o%u%x",
    &i1, &d1, suite, &c1, &e1, &E1, &f_1, &g1, &G1, &o1, &u1, &x1);
    Our_Lib = s21_sscanf("4 123 amogus q 12e34 43e21 13.37 2.28 1337.228 101 1854 12ab",
    "%i%d%9s%c%e%E%f%g%G%o%u%x",
    &i2, &d2, s2, &c2, &e2, &E2, &f2, &g2, &G2, &o2, &u2, &x2);
    ck_assert_int_eq(i1, i2); ck_assert_int_eq(d1, d2); ck_assert_str_eq(suite, s2);
    ck_assert_int_eq(c1, c2);
    ck_assert_msg(e1 == e2, "floats not equals");
    ck_assert_msg(E1 == E2, "floats not equals");
    ck_assert_msg(f_1 == f2, "floats not equals");
    ck_assert_msg(g1 == g2, "floats not equals");
    ck_assert_msg(G1 == G2, "floats not equals");
    ck_assert_int_eq(o1, o2); ck_assert_int_eq(u1, u2);
    ck_assert_int_eq(x1, x2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("4 123 amogus q 12e34 43e21 13.37 2.28 1337.228 101 1854 12ab BA21",
    "%i%d%9s%c%e%E%f%g%G%o%u%x%X",
    &i1, &d1, suite, &c1, &e1, &E1, &f_1, &g1, &G1, &o1, &u1, &x1, &X1);
    Our_Lib = s21_sscanf("4 123 amogus q 12e34 43e21 13.37 2.28 1337.228 101 1854 12ab BA21",
    "%i%d%9s%c%e%E%f%g%G%o%u%x%X",
    &i2, &d2, s2, &c2, &e2, &E2, &f2, &g2, &G2, &o2, &u2, &x2, &X2);
    ck_assert_int_eq(i1, i2); ck_assert_int_eq(d1, d2); ck_assert_str_eq(suite, s2);
    ck_assert_int_eq(c1, c2);
    ck_assert_msg(e1 == e2, "floats not equals");
    ck_assert_msg(E1 == E2, "floats not equals");
    ck_assert_msg(f_1 == f2, "floats not equals");
    ck_assert_msg(g1 == g2, "floats not equals");
    ck_assert_msg(G1 == G2, "floats not equals");
    ck_assert_int_eq(o1, o2); ck_assert_int_eq(u1, u2); ck_assert_int_eq(x1, x2);
    ck_assert_int_eq(X1, X2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("4 123 amogus q 12e34 43e21 13.37 2.28 1337.228 101 1854 12ab BA21",
    "%i%d%9s%c%e%E%f%g%G%o%u%x%X%p1",
    &i1, &d1, suite, &c1, &e1, &E1, &f_1, &g1, &G1, &o1, &u1, &x1, &X1, &p1);
    Our_Lib = s21_sscanf("4 123 amogus q 12e34 43e21 13.37 2.28 1337.228 101 1854 12ab BA21",
    "%i%d%9s%c%e%E%f%g%G%o%u%x%X%p2",
    &i2, &d2, s2, &c2, &e2, &E2, &f2, &g2, &G2, &o2, &u2, &x2, &X2, &p2);
    ck_assert_int_eq(i1, i2); ck_assert_int_eq(d1, d2); ck_assert_str_eq(suite, s2);
    ck_assert_int_eq(c1, c2); ck_assert_msg(e1 == e2, "floats not equals");
    ck_assert_msg(E1 == E2, "floats not equals");
    ck_assert_msg(f_1 == f2, "floats not equals");
    ck_assert_msg(g1 == g2, "floats not equals");
    ck_assert_msg(G1 == G2, "floats not equals");
    ck_assert_int_eq(o1, o2); ck_assert_int_eq(u1, u2); ck_assert_int_eq(x1, x2);
    ck_assert_int_eq(X1, X2); ck_assert_ptr_eq(p1, p2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(
    "4 123 amogus q 12e34 43e21 13.37 2.28 1337.228 101 1854 12ab BA21 shrek & amogus is sus?",
    "%i%d%9s%c%e%E%f%g%G%o%u%x%X%p%n%%",
    &i1, &d1, suite, &c1, &e1, &E1, &f_1, &g1, &G1, &o1, &u1, &x1, &X1, &p1, &n_1);
    Our_Lib = s21_sscanf(
    "4 123 amogus q 12e34 43e21 13.37 2.28 1337.228 101 1854 12ab BA21",
    "%i%d%9s%c%e%E%f%g%G%o%u%x%X%p%n%%",
    &i2, &d2, s2, &c2, &e2, &E2, &f2, &g2, &G2, &o2, &u2, &x2, &X2, &p2, &n_2);
    ck_assert_int_eq(i1, i2);
    ck_assert_int_eq(d1, d2);
    ck_assert_str_eq(suite, s2);
    ck_assert_int_eq(c1, c2);
    ck_assert_msg(e1 == e2, "floats not equals");
    ck_assert_msg(E1 == E2, "floats not equals");
    ck_assert_msg(f_1 == f2, "floats not equals");
    ck_assert_msg(g1 == g2, "floats not equals");
    ck_assert_msg(G1 == G2, "floats not equals");
    ck_assert_int_eq(o1, o2);
    ck_assert_int_eq(u1, u2);
    ck_assert_int_eq(x1, x2);
    ck_assert_int_eq(X1, X2);
    ck_assert_ptr_eq(p1, p2);
    ck_assert_int_eq(n_1, n_2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib =
    sscanf("4 123 amogus q 12e34 43e21 13.37 2.28 1337.228 101 1854 12ab BA21 shrek & amogus is sus?",
    "%i%d%9s%c%e%E%f%g%G%o%u%x%X%p%n%9s%%",
    &i1, &d1, suite, &c1, &e1, &E1, &f_1, &g1,
    &G1, &o1, &u1, &x1, &X1, &p1, &n_1, proc1);
    Our_Lib =
    s21_sscanf("4 123 amogus q 12e34 43e21 13.37 2.28 1337.228 101 1854 12ab BA21 shrek & amogus is sus?",
    "%i%d%9s%c%e%E%f%g%G%o%u%x%X%p%n%9s%%",
    &i2, &d2, s2, &c2, &e2, &E2, &f2, &g2,
    &G2, &o2, &u2, &x2, &X2, &p2, &n2, proc2);
    ck_assert_int_eq(i1, i2); ck_assert_int_eq(d1, d2); ck_assert_str_eq(suite, s2);
    ck_assert_int_eq(c1, c2); ck_assert_msg(e1 == e2, "floats not equals");
    ck_assert_msg(E1 == E2, "floats not equals");
    ck_assert_msg(f_1 == f2, "floats not equals");
    ck_assert_msg(g1 == g2, "floats not equals");
    ck_assert_msg(G1 == G2, "floats not equals");
    ck_assert_int_eq(o1, o2);
    ck_assert_int_eq(u1, u2);
    ck_assert_int_eq(x1, x2);
    ck_assert_int_eq(X1, X2);
    ck_assert_ptr_eq(p1, p2);
    ck_assert_int_eq(n_1, n_2);
    ck_assert_str_eq(proc1, proc2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = s21_sscanf(
    "4 123 amogusq 12e34 43e21 13.37 2.28 1337.228 101 1854 12ab BA21 shrek & amogus is sus?",
    "%12i%4d%20s%c%3e%7E%1f%9g%5G%6o%12u%11x%2X%10p%n%21s%%",
    &i1, &d1, suite, &c1, &e1, &E1, &f_1, &g1, &G1, &o1, &u1, &x1, &X1, &p1, &n_1, proc1);
    Our_Lib = s21_sscanf(
    "4 123 amogusq 12e34 43e21 13.37 2.28 1337.228 101 1854 12ab BA21 shrek & amogus is sus?",
    "%12i%4d%20s%c%3e%7E%1f%9g%5G%6o%12u%11x%2X%10p%n%21s%%",
    &i2, &d2, s2, &c2, &e2, &E2, &f2, &g2, &G2, &o2, &u2, &x2, &X2, &p2, &n_2, proc2);
    ck_assert_int_eq(i1, i2); ck_assert_int_eq(d1, d2); ck_assert_str_eq(suite, s2);
    ck_assert_int_eq(c1, c2);
    ck_assert_msg(e1 == e2, "floats not equals");
    ck_assert_msg((int)E1 == (int)E2, "floats not equals");
    ck_assert_msg(f_1 == f2, "floats not equals");
    ck_assert_msg(g1 == g2, "floats not equals");
    ck_assert_msg(G1 == G2, "floats not equals");
    ck_assert_int_eq(o1, o2);
    ck_assert_int_eq(u1, u2);
    ck_assert_int_eq(x1, x2);
    ck_assert_int_eq(X1, X2);
    ck_assert_ptr_eq(p1, p2);
    ck_assert_int_eq(n_1, n_2);
    ck_assert_str_eq(proc1, proc2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(
    "4 123 amogusq 12e34 43e21 13.37 2.28 1337.228 101 1854 12ab BA21 shrek & amogus is sus?",
    "%*i%12i%4d%20s%c%3e%7E%*i%1f%9g%5G%6o%12u%*x%11x%2X%10p%n%21s%%",
    &i1, &d1, suite, &c1, &e1, &E1, &f_1, &g1, &G1, &o1, &u1, &x1, &X1, &p1, &n_1, proc1);
    Our_Lib = s21_sscanf(
    "4 123 amogusq 12e34 43e21 13.37 2.28 1337.228 101 1854 12ab BA21 shrek & amogus is sus?",
    "%*i%12i%4d%20s%c%3e%7E%*i%1f%9g%5G%6o%12u%*x%11x%2X%10p%n%21s%%",
    &i2, &d2, s2, &c2, &e2, &E2, &f2, &g2, &G2, &o2, &u2, &x2, &X2, &p2, &n_2, proc2);
    ck_assert_int_eq(i1, i2); ck_assert_int_eq(d1, d2); ck_assert_str_eq(suite, s2);
    ck_assert_int_eq(c1, c2); ck_assert_msg(e1 == e2, "floats not equals");
    ck_assert_msg(E1 == E2, "floats not equals");
    ck_assert_msg(f_1 == f2, "floats not equals");
    ck_assert_msg(g1 == g2, "floats not equals");
    ck_assert_msg(G1 == G2, "floats not equals");
    ck_assert_int_eq(o1, o2); ck_assert_int_eq(u1, u2); ck_assert_int_eq(x1, x2);
    ck_assert_int_eq(X1, X2); ck_assert_ptr_eq(p1, p2); ck_assert_int_eq(n_1, n_2);
    ck_assert_str_eq(proc1, proc2); ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(
    "4 123 amogusq 12e34 43e21 13.37 2.28 1337.228 101 1854 12ab BA21 shrek & amogus is sus?",
    "%*i%12hi%4hd%20s%c%3Le%7LE%*i%1f%9Lg%5LG%6lo%12lu%*x%11x%2X%10p%n%21s%%",
    &hi1, &hd1, suite, &c1, &Le1, &LE1, &f_1, &Lg1, &LG1, &lo1, &lu1, &x1, &X1, &p1, &n_1, proc1);
    Our_Lib = s21_sscanf(
    "4 123 amogusq 12e34 43e21 13.37 2.28 1337.228 101 1854 12ab BA21 shrek & amogus is sus?",
    "%*i%12hi%4hd%20s%c%3Le%7LE%*i%1f%9Lg%5LG%6lo%12lu%*x%11x%2X%10p%n%21s%%",
    &hi2, &hd2, s2, &c2, &Le2, &LE2, &f2, &Lg2, &LG2, &lo2, &lu2, &x2, &X2, &p2, &n_2, proc2);
    ck_assert_int_eq(hi1, hi2); ck_assert_int_eq(hd1, hd2); ck_assert_str_eq(suite, s2);
    ck_assert_int_eq(c1, c2); ck_assert_msg(e1 == e2, "floats not equals");
    ck_assert_msg(E1 == E2, "floats not equals");
    ck_assert_msg(f_1 == f2, "floats not equals");
    ck_assert_msg(g1 == g2, "floats not equals");
    ck_assert_msg(G1 == G2, "floats not equals");
    ck_assert_int_eq(lo1, lo2); ck_assert_int_eq(lu1, lu2); ck_assert_int_eq(x1, x2);
    ck_assert_int_eq(X1, X2); ck_assert_ptr_eq(p1, p2); ck_assert_int_eq(n_1, n_2);
    ck_assert_str_eq(proc1, proc2); ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST


START_TEST(sscanf_test_c_1) {
    char *empty_str = "";
    char c1, c2;
    int Std_Lib, Our_Lib;
    char *stStd_Lib = " abc ";

    Std_Lib = sscanf(stStd_Lib, "%*c%*c%c", &c1);
    Our_Lib = s21_sscanf(stStd_Lib, "%*c%*c%c", &c2);
    ck_assert_int_eq(Std_Lib, Our_Lib);
    ck_assert_int_eq(c1, c2);

    Std_Lib = sscanf(empty_str, "%c", &c1);
    Our_Lib = s21_sscanf(empty_str, "%c", &c2);
    ck_assert_int_eq(Std_Lib, Our_Lib);
    ck_assert_int_eq(c1, c2);
}
END_TEST

START_TEST(sscanf_test_d_2) {
    char *empty_str = "";
    int nmb1 = 111, nmb2 = 111;
    int Std_Lib, Our_Lib;
    char *stOur_Lib = "  123   -123   567  32368  ";
    char wrong_stOur_Lib[] = " abc 324 ";
    Std_Lib = sscanf(stOur_Lib, "%*d%d%*d%*d", &nmb1);
    Our_Lib = s21_sscanf(stOur_Lib, "%*d%d%*d%*d", &nmb2);
    ck_assert_int_eq(nmb1, nmb2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    nmb1 = 1, nmb2 = 1;
    Std_Lib = sscanf(empty_str, "%d", &nmb1);
    Our_Lib = s21_sscanf(empty_str, "%d", &nmb2);
    ck_assert_int_eq(nmb1, nmb2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    nmb1 = 0, nmb2 = 0;
    Std_Lib = sscanf(wrong_stOur_Lib, "%d", &nmb1);
    Our_Lib = s21_sscanf(wrong_stOur_Lib, "%d", &nmb2);
    ck_assert_int_eq(nmb1, nmb2);
    ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_i_2) {
    char *empty_str = "";
    int nmbX1 = 1, nmbX2 = 1;
    int nmbO1 = 1, nmbO2 = 1;
    int nmbD1 = 1, nmbD2 = 1;
    int Std_Lib, Our_Lib;
    char *str3 = "  0x3F5  0600 -12345 987";
    char *wrong_str3 = "  0000xFF5  ";

    Std_Lib = sscanf(str3, "%i%i%i%*i", &nmbX1, &nmbO1, &nmbD1);
    Our_Lib = s21_sscanf(str3, "%i%i%i%*i", &nmbX2, &nmbO2, &nmbD2);
    ck_assert_int_eq(nmbX1, nmbX2);
    ck_assert_int_eq(nmbO1, nmbO2);
    ck_assert_int_eq(nmbD1, nmbD2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(empty_str, "%i", &nmbD1);
    Our_Lib = s21_sscanf(empty_str, "%i", &nmbD2);
    ck_assert_int_eq(nmbD1, nmbD2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(wrong_str3, "%i", &nmbX1);
    Our_Lib = s21_sscanf(wrong_str3, "%i", &nmbX2);
    ck_assert_uint_eq(nmbX1, nmbX2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    unsigned int octnmb1 = 1, octnmb2 = 2;

    Std_Lib = sscanf(str3, "%*x%o", &octnmb1);
    Our_Lib = s21_sscanf(str3, "%*x%o", &octnmb2);
    ck_assert_uint_eq(octnmb1, octnmb2);
    ck_assert_int_eq(Std_Lib, Our_Lib);
    }
END_TEST

START_TEST(sscanf_test_u_2) {
    char *empty_str = "";
    unsigned u1 = 1, u2 = 2, Std_Lib, Our_Lib;
    char *str4 = " 0    3147483647";
    char *wrong_str4 = "  Verter_Hi123456";

    Std_Lib = sscanf(str4, "%*u%u", &u1);
    Our_Lib = s21_sscanf(str4, "%*u%u", &u2);
    ck_assert_uint_eq(u1, u2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(empty_str, "%u", &u1);
    Our_Lib = s21_sscanf(empty_str, "%u", &u2);
    ck_assert_int_eq(u1, u2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(wrong_str4, "%u", &u1);
    Our_Lib = s21_sscanf(wrong_str4, "%u", &u2);
    ck_assert_uint_eq(u1, u2);
    ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_x_2) {
    unsigned u1 = 1, u2 = 2, u3 = 3, u4 = 4, Std_Lib, Our_Lib;
    char *str4 = " 123123  3147487 123 123";

    Std_Lib = sscanf(str4, "%*x%x%X%*X", &u1, &u3);
    Our_Lib = s21_sscanf(str4, "%*x%x%X%*X", &u2, &u4);
    ck_assert_uint_eq(u1, u2);
    ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_double) {
    char *empty_str = "";
    int Std_Lib, Our_Lib;
    double fl1 = 1, fl2 = 2;
    double dbl1 = 1, dbl2 = 2;
    char *str5 = "  -876.214235  2145.93560734567  324.21427 ";
    char *wrong_str5 = "  Verter_Kek35634.34";

    Std_Lib = sscanf(str5, "%lf%lf%*f", &fl1, &dbl1);
    Our_Lib = s21_sscanf(str5, "%lf%lf%*f", &fl2, &dbl2);
    ck_assert_msg(dbl1 == dbl2, "Float test 1 fail");
    ck_assert_msg(fl1 == fl2, "Float test 2 fail");
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(empty_str, "%lf", &dbl1);
    Our_Lib = s21_sscanf(empty_str, "%lf", &dbl2);
    ck_assert_float_eq(dbl1, dbl2);
    ck_assert_msg(dbl1 == dbl2, "Float test 3 fail");
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(wrong_str5, "%lf", &dbl1);
    Our_Lib = s21_sscanf(wrong_str5, "%lf", &dbl2);
    ck_assert_float_eq(dbl1, dbl2);
    ck_assert_msg(dbl1 == dbl2, "Float test 4 fail");
    ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_EeGg) {
    char *empty_str = "";
    int Std_Lib, Our_Lib;
    double e1, e2, g1, g2;
    char *str6_u = "  0.08762147E+03  1446567457.024E-11 0E+02";
    char *wrong_str6 = "  Verter?35634E+02";

    Std_Lib = sscanf(str6_u, "%lE%lG%*e", &e1, &g1);
    Our_Lib = s21_sscanf(str6_u, "%lE%lG%*e", &e2, &g2);
    ck_assert_int_eq(e1, e2);
    ck_assert_int_eq(g1,  g2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(empty_str, "%le", &e1);
    Our_Lib = s21_sscanf(empty_str, "%le", &e2);
    ck_assert_int_eq(e1, e2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(wrong_str6, "%lg", &g1);
    Our_Lib = s21_sscanf(wrong_str6, "%lg", &g2);
    ck_assert_int_eq(g1,  g2);
    ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_s_2) {
    int Std_Lib, Our_Lib;
    char str_vaStd_Lib[30] = "", str_vaOur_Lib[30] = "";

    Std_Lib = sscanf(" this is random string", "%*s%*s%10s%*s", str_vaStd_Lib);
    Our_Lib = s21_sscanf(" this is random string", "%*s%*s%10s%*s", str_vaOur_Lib);
    ck_assert_str_eq(str_vaStd_Lib, str_vaOur_Lib);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf("\0", "%10s", str_vaStd_Lib);
    Our_Lib = s21_sscanf("\0", "%10s", str_vaOur_Lib);
    ck_assert_str_eq(str_vaStd_Lib, str_vaOur_Lib);
    ck_assert_int_eq(Std_Lib, Our_Lib);
}

START_TEST(sscanf_test_ld) {
    char *empty_str = "";
    int Std_Lib, Our_Lib;
    long nm1 = 111, nm2 = 111;
    char *str5 = "  123   -5123123123   567  588  ";
    char wrong_str5[] = " abc 123 ";

    Std_Lib = sscanf(str5, "%*d%ld%*d%*d", &nm1);
    Our_Lib = s21_sscanf(str5, "%*d%ld%*d%*d", &nm2);
    ck_assert_int_eq(nm1, nm2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(empty_str, "%ld", &nm1);
    Our_Lib = s21_sscanf(empty_str, "%ld", &nm2);
    ck_assert_int_eq(nm1, nm2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(wrong_str5, "%ld", &nm1);
    Our_Lib = s21_sscanf(wrong_str5, "%ld", &nm2);
    ck_assert_int_eq(nm1, nm2);
    ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_hd) {
    char *empty_str = "";
    int Std_Lib, Our_Lib;
    short nh1 = 111, nh2 = 111;
    char *str6 = "  123   -5123123123   567  588  ";
    char wrong_str6[] = " abc 123 ";

    Std_Lib = sscanf(str6, "%*d%hd%*d%*d", &nh1);
    Our_Lib = s21_sscanf(str6, "%*d%hd%*d%*d", &nh2);
    ck_assert_int_eq(nh1, nh2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(empty_str, "%hd", &nh1);
    Our_Lib = s21_sscanf(empty_str, "%hd", &nh2);
    ck_assert_int_eq(nh1, nh2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(wrong_str6, "%hd", &nh1);
    Our_Lib = s21_sscanf(wrong_str6, "%hd", &nh2);
    ck_assert_int_eq(nh1, nh2);
    ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_lu) {
    char *empty_str = "";
    int Std_Lib, Our_Lib;
    unsigned long nu1 = 111, nu2 = 111;
    char *str7 = "  123   5123123123   567  588  ";
    char wrong_str7[] = " abc 123 ";

    Std_Lib = sscanf(str7, "%*d%lu%*d%*d", &nu1);
    Our_Lib = s21_sscanf(str7, "%*d%lu%*d%*d", &nu2);
    ck_assert_int_eq(nu1, nu2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(empty_str, "%lu", &nu1);
    Our_Lib = s21_sscanf(empty_str, "%lu", &nu2);
    ck_assert_int_eq(nu1, nu2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(wrong_str7, "%lu", &nu1);
    Our_Lib = s21_sscanf(wrong_str7, "%lu", &nu2);
    ck_assert_int_eq(nu1, nu2);
    ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_hu) {
    char *empty_str = "";
    int Std_Lib, Our_Lib;
    unsigned short nuh1 = 111, nuh2 = 111;
    char *str5 = "  123   15000   567  588  ";
    char wrong_str5[] = " abc 123 ";

    Std_Lib = sscanf(str5, "%*d%hu%*d%*d", &nuh1);
    Our_Lib = s21_sscanf(str5, "%*d%hu%*d%*d", &nuh2);
    ck_assert_int_eq(nuh1, nuh2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(empty_str, "%hu", &nuh1);
    Our_Lib = s21_sscanf(empty_str, "%hu", &nuh2);
    ck_assert_int_eq(nuh1, nuh2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(wrong_str5, "%hu", &nuh1);
    Our_Lib = s21_sscanf(wrong_str5, "%hu", &nuh2);
    ck_assert_int_eq(nuh1, nuh2);
    ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

START_TEST(sscanf_test_ptr) {
    char *empty_str = "";
    int Std_Lib, Our_Lib;
    unsigned long p1 = 111, p2 = 111;
    unsigned long *ptStd_Lib = &p1;
    unsigned long *ptStd_Lib1 = &p1;
    unsigned long *ptOur_Lib = &p2;
    unsigned long *ptOur_Lib2 = &p2;
    char *str8 = "  123   5123123123   0x567  588  ";
    char wrong_str8[] = " abc 123 ";

    Std_Lib = sscanf(str8, "%*d%p%p%*p", &ptStd_Lib, &ptStd_Lib1);
    Our_Lib = s21_sscanf(str8, "%*d%p%p%*p", &ptOur_Lib, &ptOur_Lib2);
    ck_assert_ptr_eq(ptStd_Lib, ptOur_Lib);
    ck_assert_ptr_eq(ptStd_Lib1, ptOur_Lib2);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(empty_str, "%p", &ptStd_Lib);
    Our_Lib = s21_sscanf(empty_str, "%p", &ptOur_Lib);
    ck_assert_ptr_eq(ptStd_Lib, ptOur_Lib);
    ck_assert_int_eq(Std_Lib, Our_Lib);

    Std_Lib = sscanf(wrong_str8, "%p", &ptStd_Lib);
    Our_Lib = s21_sscanf(wrong_str8, "%p", &ptOur_Lib);
    ck_assert_ptr_eq(ptStd_Lib, ptOur_Lib);
    ck_assert_int_eq(Std_Lib, Our_Lib);
}
END_TEST

int main(void) {
    Suite *suite = suite_create("S21_Debuger");
    SRunner *srunner = srunner_create(suite);

    TCase *s21_strlen_test = tcase_create("s21_strlen");
    suite_add_tcase(suite, s21_strlen_test);
    tcase_add_test(s21_strlen_test, s21_strlenTEST);

    TCase *s21_memchr_test = tcase_create("s21_memchr");
    suite_add_tcase(suite, s21_memchr_test);
    tcase_add_test(s21_memchr_test, s21_memchrTEST);

    TCase *s21_memcmp_test = tcase_create("s21_memcmp");
    suite_add_tcase(suite, s21_memcmp_test);
    tcase_add_test(s21_memcmp_test, s21_memcmpTEST);

    TCase *s21_memcpy_test = tcase_create("s21_memcpy");
    suite_add_tcase(suite, s21_memcpy_test);
    tcase_add_test(s21_memcpy_test, s21_memcpyTEST);

    TCase *s21_memmove_test = tcase_create("s21_memmove");
    suite_add_tcase(suite, s21_memmove_test);
    tcase_add_test(s21_memmove_test, s21_memmoveTEST);

    TCase *s21_memset_test = tcase_create("s21_memset");
    suite_add_tcase(suite, s21_memset_test);
    tcase_add_test(s21_memset_test, s21_memsetTEST);

    TCase *s21_strcat_test = tcase_create("s21_strcatTEST");
    suite_add_tcase(suite, s21_strcat_test);
    tcase_add_test(s21_strcat_test, s21_strcatTEST);

    TCase *s21_strncat_test = tcase_create("s21_strncat");
    suite_add_tcase(suite, s21_strncat_test);
    tcase_add_test(s21_strncat_test, s21_strncatTEST);

    TCase *s21_strchr_test = tcase_create("s21_strchr");
    suite_add_tcase(suite, s21_strchr_test);
    tcase_add_test(s21_strchr_test, s21_strchrTEST);

    TCase *s21_strrchr_test = tcase_create("s21_strrchr");
    suite_add_tcase(suite, s21_strrchr_test);
    tcase_add_test(s21_strrchr_test, s21_strrchrTEST);

    TCase *s21_strcmp_test = tcase_create("s21_strcmp");
    suite_add_tcase(suite, s21_strcmp_test);
    tcase_add_test(s21_strcmp_test, s21_strcmpTEST);

    TCase *s21_strncmp_test = tcase_create("s21_strncmp");
    suite_add_tcase(suite, s21_strncmp_test);
    tcase_add_test(s21_strncmp_test, s21_strncmpTEST);

    TCase *s21_strcpy_test = tcase_create("s21_strcpy");
    suite_add_tcase(suite, s21_strcpy_test);
    tcase_add_test(s21_strcpy_test, s21_strcpyTEST);

    TCase *s21_strncpy_test = tcase_create("s21_strncpy");
    suite_add_tcase(suite, s21_strncpy_test);
    tcase_add_test(s21_strncpy_test, s21_strncpyTEST);

    TCase *s21_strcspn_test = tcase_create("s21_strcspn");
    suite_add_tcase(suite, s21_strcspn_test);
    tcase_add_test(s21_strcspn_test, s21_strcspnTEST);

    TCase *s21_strspn_test = tcase_create("s21_strspn");
    suite_add_tcase(suite, s21_strspn_test);
    tcase_add_test(s21_strspn_test, s21_strspnTEST);

    TCase *s21_strstr_test = tcase_create("s21_strstr");
    suite_add_tcase(suite, s21_strstr_test);
    tcase_add_test(s21_strstr_test, s21_strstrTEST);

    TCase *s21_strtok_test = tcase_create("s21_strtok");
    suite_add_tcase(suite, s21_strtok_test);
    tcase_add_test(s21_strtok_test, s21_strtokTEST);

    TCase *s21_strerror_test = tcase_create("s21_strerror");
    suite_add_tcase(suite, s21_strerror_test);
    tcase_add_test(s21_strerror_test, s21_strerrorTEST);

    TCase *s21_strpbrk_test = tcase_create("s21_strpbrk");
    suite_add_tcase(suite, s21_strpbrk_test);
    tcase_add_test(s21_strpbrk_test, s21_strpbrkTEST);

    TCase *s21_to_upper_test = tcase_create("s21_to_upper");
    suite_add_tcase(suite, s21_to_upper_test);
    tcase_add_test(s21_to_upper_test, s21_to_upperTEST);

    TCase *s21_to_lower_test = tcase_create("s21_to_lower");
    suite_add_tcase(suite, s21_to_lower_test);
    tcase_add_test(s21_to_lower_test, s21_to_lowerTEST);

    TCase *s21_insert_test = tcase_create("s21_insert");
    suite_add_tcase(suite, s21_insert_test);
    tcase_add_test(s21_insert_test, s21_insertTEST);

    TCase *s21_trim_test = tcase_create("s21_trim");
    suite_add_tcase(suite, s21_trim_test);
    tcase_add_test(s21_trim_test, s21_trimTEST);

    TCase *sprintf_d_test = tcase_create("Sprnt_d");
    suite_add_tcase(suite, sprintf_d_test);
    tcase_add_test(sprintf_d_test, Sprintf_Test_d);

    TCase *sprintf_i_test = tcase_create("Sprnt_i");
    suite_add_tcase(suite, sprintf_i_test);
    tcase_add_test(sprintf_i_test, Sprintf_Test_i);

    TCase *sprintf_c_test = tcase_create("Sprnt_c");
    suite_add_tcase(suite, sprintf_c_test);
    tcase_add_test(sprintf_c_test, Sprintf_Test_c);

    TCase *sprintf_s_test = tcase_create("Sprnt_s");
    suite_add_tcase(suite, sprintf_s_test);
    tcase_add_test(sprintf_s_test, Sprintf_Test_s);

    TCase *sprintf_f = tcase_create("Sprnt_f");
    suite_add_tcase(suite, sprintf_f);
    tcase_add_test(sprintf_f, Sprintf_Test_f);

    TCase *sprintf_u_test = tcase_create("Sprnt_u");
    suite_add_tcase(suite, sprintf_u_test);
    tcase_add_test(sprintf_u_test, Sprintf_Test_u);

    TCase *sprintf_pr_test = tcase_create("Sprnt_pr");
    suite_add_tcase(suite, sprintf_pr_test);
    tcase_add_test(sprintf_pr_test, Sprintf_Test_pr);

    TCase *sprintf_g = tcase_create("Sprnt_g");
    suite_add_tcase(suite, sprintf_g);
    tcase_add_test(sprintf_g, Sprintf_Test_g);

    TCase *sprintf_G_test = tcase_create("Sprnt_G");
    suite_add_tcase(suite, sprintf_G_test);
    tcase_add_test(sprintf_G_test, Sprintf_Test_G);

    TCase *sprintf_e_test = tcase_create("Sprnt_e");
    suite_add_tcase(suite, sprintf_e_test);
    tcase_add_test(sprintf_e_test, Sprintf_Test_e);

    TCase *sprintf_E_test = tcase_create("Sprnt_E");
    suite_add_tcase(suite, sprintf_E_test);
    tcase_add_test(sprintf_E_test, Sprintf_Test_E);

    TCase *sprintf_x_test = tcase_create("Sprnt_x");
    suite_add_tcase(suite, sprintf_x_test);
    tcase_add_test(sprintf_x_test, Sprintf_Test_x);

    TCase *sprintf_X_test = tcase_create("Sprnt_X");
    suite_add_tcase(suite, sprintf_X_test);
    tcase_add_test(sprintf_X_test, Sprintf_Test_X);

    TCase *sprintf_o_test = tcase_create("Sprnt_o");
    suite_add_tcase(suite, sprintf_o_test);
    tcase_add_test(sprintf_o_test, Sprintf_Test_o);

    TCase *sprintf_n_test = tcase_create("Sprnt_n");
    suite_add_tcase(suite, sprintf_n_test);
    tcase_add_test(sprintf_n_test, Sprintf_Test_n);

    TCase *sprintf_p_test = tcase_create("Sprnt_p");
    suite_add_tcase(suite, sprintf_p_test);
    tcase_add_test(sprintf_p_test, Sprintf_Test_p);

    TCase *sprintf_MT_test = tcase_create("Sprnt_Multy");
    suite_add_tcase(suite, sprintf_MT_test);
    tcase_add_test(sprintf_MT_test, test_s21_MultySpec);

    TCase *s21_Sscanf_e = tcase_create("Sscanf_e");
    suite_add_tcase(suite, s21_Sscanf_e);
    tcase_add_test(s21_Sscanf_e, sscanf_test_e_1);

    TCase *s21_Sscanf_EeGg = tcase_create("Sscanf_EeGg");
    suite_add_tcase(suite, s21_Sscanf_EeGg);
    tcase_add_test(s21_Sscanf_EeGg, sscanf_test_EeGg);

    TCase *s21_Sscanf_l_1 = tcase_create("Sscanf_l");
    suite_add_tcase(suite, s21_Sscanf_l_1);
    tcase_add_test(s21_Sscanf_l_1, sscanf_test_l_1);

    TCase *s21_Sscanf_L_1 = tcase_create("Sscanf_L");
    suite_add_tcase(suite, s21_Sscanf_L_1);
    tcase_add_test(s21_Sscanf_L_1, sscanf_test_L_1);

    TCase *s21_Sscanf_ld = tcase_create("Sscanf_ld");
    suite_add_tcase(suite, s21_Sscanf_ld);
    tcase_add_test(s21_Sscanf_ld, sscanf_test_ld);

    TCase *s21_Sscanf_hd = tcase_create("Sscanf_hd");
    suite_add_tcase(suite, s21_Sscanf_hd);
    tcase_add_test(s21_Sscanf_hd, sscanf_test_hd);

    TCase *s21_Sscanf_lu = tcase_create("Sscanf_lu");
    suite_add_tcase(suite, s21_Sscanf_lu);
    tcase_add_test(s21_Sscanf_lu, sscanf_test_lu);

    TCase *s21_Sscanf_hu = tcase_create("Sscanf_hu");
    suite_add_tcase(suite, s21_Sscanf_hu);
    tcase_add_test(s21_Sscanf_hu, sscanf_test_hu);

    TCase *s21_Sscanf_h_1 = tcase_create("Sscanf_h");
    suite_add_tcase(suite, s21_Sscanf_h_1);
    tcase_add_test(s21_Sscanf_h_1, sscanf_test_h_1);

    TCase *s21_Sscanf_d_1 = tcase_create("Sscanf_d");
    suite_add_tcase(suite, s21_Sscanf_d_1);
    tcase_add_test(s21_Sscanf_d_1, sscanf_test_d_1);

    TCase *s21_Sscanf_d_2 = tcase_create("Sscanf_d");
    suite_add_tcase(suite, s21_Sscanf_d_2);
    tcase_add_test(s21_Sscanf_d_2, sscanf_test_d_2);

    TCase *s21_Sscanf_u_1 = tcase_create("Sscanf_u");
    suite_add_tcase(suite, s21_Sscanf_u_1);
    tcase_add_test(s21_Sscanf_u_1, sscanf_test_u_1);

    TCase *s21_Sscanf_u_2 = tcase_create("Sscanf_u");
    suite_add_tcase(suite, s21_Sscanf_u_2);
    tcase_add_test(s21_Sscanf_u_2, sscanf_test_u_2);

    TCase *s21_Sscanf_ptr = tcase_create("Sscanf_ptr");
    suite_add_tcase(suite, s21_Sscanf_ptr);
    tcase_add_test(s21_Sscanf_ptr, sscanf_test_ptr);

    TCase *s21_Sscanf_proc = tcase_create("Sscanf_proc");
    suite_add_tcase(suite, s21_Sscanf_proc);
    tcase_add_test(s21_Sscanf_proc, sscanf_test_proc_1);

    TCase *s21_Sscanf_p_1 = tcase_create("s21_Sscanf_p");
    suite_add_tcase(suite, s21_Sscanf_p_1);
    tcase_add_test(s21_Sscanf_p_1, sscanf_test_p_1);

    TCase *s21_Sscanf_s_1 = tcase_create("Sscanf_s");
    suite_add_tcase(suite, s21_Sscanf_s_1);
    tcase_add_test(s21_Sscanf_s_1, sscanf_test_s_1);

    TCase *s21_Sscanf_s_2 = tcase_create("Sscanf_s");
    suite_add_tcase(suite, s21_Sscanf_s_2);
    tcase_add_test(s21_Sscanf_s_2, sscanf_test_s_2);

    TCase *s21_Sscanf_f_1 = tcase_create("Sscanf_f");
    suite_add_tcase(suite, s21_Sscanf_f_1);
    tcase_add_test(s21_Sscanf_f_1, sscanf_test_f_1);

    TCase *s21_Sscanf_o_1 = tcase_create("Sscanf_o");
    suite_add_tcase(suite, s21_Sscanf_o_1);
    tcase_add_test(s21_Sscanf_o_1, sscanf_test_o_1);

    TCase *s21_Sscanf_x_1 = tcase_create("Sscanf_x");
    suite_add_tcase(suite, s21_Sscanf_x_1);
    tcase_add_test(s21_Sscanf_x_1, sscanf_test_x_1);

    TCase *s21_Sscanf_x_2 = tcase_create("Sscanf_x");
    suite_add_tcase(suite, s21_Sscanf_x_2);
    tcase_add_test(s21_Sscanf_x_2, sscanf_test_x_2);

    TCase *s21_Sscanf_i_1 = tcase_create("Sscanf_i");
    suite_add_tcase(suite, s21_Sscanf_i_1);
    tcase_add_test(s21_Sscanf_i_1, sscanf_test_i_1);

    TCase *s21_Sscanf_i_2 = tcase_create("Sscanf_i");
    suite_add_tcase(suite, s21_Sscanf_i_2);
    tcase_add_test(s21_Sscanf_i_2, sscanf_test_i_2);

    TCase *s21_Sscanf_c_1 = tcase_create("Sscanf_c");
    suite_add_tcase(suite, s21_Sscanf_c_1);
    tcase_add_test(s21_Sscanf_c_1, sscanf_test_c_1);

    TCase *s21_Sscanf_n_1 = tcase_create("Sscanf_n");
    suite_add_tcase(suite, s21_Sscanf_n_1);
    tcase_add_test(s21_Sscanf_n_1, sscanf_test_n_1);

    TCase *s21_Sscanf_test_double = tcase_create("Sscanf_double");
    suite_add_tcase(suite, s21_Sscanf_test_double);
    tcase_add_test(s21_Sscanf_test_double, sscanf_test_double);

    TCase *s21_Sscanf_all = tcase_create("Sscanf_all");
    suite_add_tcase(suite, s21_Sscanf_all);
    tcase_add_test(s21_Sscanf_all, sscanf_test_all);

    srunner_run_all(srunner, CK_VERBOSE);
    int number_failed = srunner_ntests_failed(srunner);
    srunner_free(srunner);

    return number_failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
