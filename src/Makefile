OBJECTS = chpt.o chpt_R.o chpt.so

all: chpt.so

clean:
	rm -f $(OBJECTS)

chpt.so: chpt.o chpt_R.c chpt.h
	R CMD SHLIB chpt_R.c chpt.o -o chpt.so -lm

chpt.o: chpt.c chpt.h
	gcc -std=gnu99 -fPIC -g -O0 -c chpt.c \
		-lpthread -Wall
