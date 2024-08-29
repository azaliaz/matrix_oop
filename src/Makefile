CC=g++
FLAGS= -Wall -Werror -Wextra -std=c++17
LIBS=-lgtest -lstdc++
OS:=$(shell uname -s)

ifeq ($(OS), Darwin)
	OPEN=open
	LEAKS=leaks --atExit --
else
	LIBS+=-lm
	OPEN=xdg-open
	LEAKS=valgrind
endif

REPORT_DIR = report

rebuild: clean all

all: test

s21_matrix_oop.a: s21_matrix_oop.cpp
	$(CC) $(FLAGS) -c s21_matrix_oop.cpp -o s21_matrix_oop.o
	ar rc s21_matrix_oop.a s21_matrix_oop.o
	ranlib s21_matrix_oop.a

test: s21_matrix_oop_test.cpp s21_matrix_oop.a
	$(CC) $(FLAGS) -c s21_matrix_oop_test.cpp -o s21_matrix_oop_test.o
	$(CC) $(FLAGS) s21_matrix_oop_test.o s21_matrix_oop.a $(LIBS) -o test
	./test

gcov_report: test
	lcov -t "gcov_report" --ignore-errors mismatch -o report.info --no-external -c -d .
	genhtml -o report report.info
	open ./report/index.html

clean:
	rm -rf *.o *.a *.gcda *.gcno *.info ./report test

clf:
	clang-format -n --style=Google *.cpp *.h

cpk:
	# cppcheck --language=c++ --check-level=exhaustive *.cpp *.h
	cppcheck --language=c++ *.cpp *.h

style: clf cpk

leaks: test
	$(LEAKS) ./test
