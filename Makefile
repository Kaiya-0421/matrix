bl: mat_ch_bl.cpp
	g++ mat_ch_bl.cpp -I /usr/local/include/eigen3 -std=c++14 -O3 -o $@
bl2: block.cpp
	g++ block.cpp -O3 -o $@
omp: omp.cpp
	g++ omp.cpp -I /usr/local/include/eigen3 -std=c++14 -fopenmp -O3 -o $@
clean:
	\rm -rf a.out fn mat bl bl2 
p: pthread.cpp
	g++ pthread.cpp -I /usr/local/include/eigen3 -std=c++14 -fopenmp -lpthread -O3 -o $@
p4: pthread4.cpp
	g++ pthread4.cpp -I /usr/local/include/eigen3 -std=c++14 -fopenmp -lpthread -O3 -o $@
p8: pthread8.cpp
	g++ pthread8.cpp -I /usr/local/include/eigen3 -std=c++14 -fopenmp -lpthread -O3 -o $@
p16: pthread16.cpp
	g++ pthread16.cpp -I /usr/local/include/eigen3 -std=c++14 -fopenmp -lpthread -O3 -o $@