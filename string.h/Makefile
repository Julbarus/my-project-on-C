WWW = -Wall -Wextra -Werror -std=c11
LINUX=-lsubunit -lrt -lpthread -lm
LIBS=-lcheck
LFLAGS=-fprofile-arcs -ftest-coverage
OS=$(shell uname -s)

all: clean test


test: s21_string.a
ifeq ($(OS), Darwin)
	gcc $(WWW) unit_test.c s21_string.a -o unit_test $(LIBS)
	make clean_obj
else
	gcc $(WWW) unit_test.c s21_string.a -o unit_test $(LIBS) $(LINUX)
	make clean_obj
endif
	./unit_test

gcov_report: clean_test
ifeq ($(OS), Darwin)
	gcc $(WWW) -fprofile-arcs -ftest-coverage unit_test.c s21_*.c -o unit_test $(LIBS)
else
	gcc $(WWW) -fprofile-arcs -ftest-coverage unit_test.c s21_*.c -o unit_test $(LIBS) $(LINUX)
endif
	./unit_test
	lcov -t "test" -o test.info -c -d .
	genhtml -o report test.info

s21_string.a: build_s21_string clean_obj

build_s21_string:
	gcc $(WWW) -c s21_*.c
	ar rcs s21_string.a s21_*.o
	ranlib s21_string.a

clean_obj:
	rm -rf *.o

clean_lib: 
	rm -rf *.a

clean_test:
	rm -rf *.gcda
	rm -rf *.gcno
	rm -rf *.info
	rm -rf test
	rm -rf report

clean: clean_lib clean_lib clean_test clean_obj
	rm -rf unit_test
