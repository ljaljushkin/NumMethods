all:
	mkdir -p bin
	icpc src/*.cpp -I include -o bin/solution -O2 -lrt -openmp

run:
	cp ../test_data/input bin/
	cd bin; ./solution input output time; cat time; echo

clean:
	rm -rvf *.o bin/* __cmake
