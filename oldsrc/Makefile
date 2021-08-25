CC = gcc
CFLAGS = -WALL -I./include
TARGET = UILC_Calculation
execute: main.o UILC_CA.o UILC_Core.o
	gcc -o exe
main.o: min.o


# gcc -o ./object/UILC_Morena.o -c -Wall ./src/UILC_Morena.c -lgsl -lm
# gcc -o ./object/UILC_General.o -c -Wall ./src/UILC_General.c -lgsl -lm
# gcc -o ./bin/test3.out -Wall ./src/test.c ./src/UILC_General.c ./src/UILC_Morena.c -lm -lgsl -lgslcblas