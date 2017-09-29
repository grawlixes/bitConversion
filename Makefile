test: floatx
	./floatx 32 8 <test1.txt

floatx:  	main.c floatx.c floatx.h
		gcc  -g -Wall -std=c99 -lm  -o floatx main.c floatx.c 

clean:
		-rm floatx
