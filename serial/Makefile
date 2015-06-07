all:
	mkdir -p bin
	icpc src/*.cpp -I include -o bin/solution -O2 -lrt

run:
	cd bin; ./solution input output time

clean:
	rm -rvf *.o bin/solution bin/output bin/time output time __cmake
