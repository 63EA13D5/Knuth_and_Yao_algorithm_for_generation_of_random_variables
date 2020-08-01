LIB_DIR=/usr/local/lib
NTL_INC=/usr/local/include/NTL
GMPL_INC=/usr/local/include

all: clean
	g++ -Wall -Wpedantic -g -O5 -I${NTL_INC} -I${GMPL_INC} ky_v4.cpp -o ky_v4 -std=c++11 -L${LIB_DIR} -lntl -lgmp -lm -lpthread

ky_v4: clean
	g++ -Wall -Wpedantic -g -O5 -I${NTL_INC} -I${GMPL_INC} ky_v4.cpp -o ky_v4 -std=c++11 -L${LIB_DIR} -lntl -lgmp -lm -lpthread

debug: clean
	g++ -Wall -Wpedantic -g -I${NTL_INC} -I${GMPL_INC} ky_v4.cpp -o ky_v4 -std=c++11 -L${LIB_DIR} -lntl -lgmp -lm -lpthread
	gdb hanmat_mt

clean:
	rm -rf ky_v4 test
