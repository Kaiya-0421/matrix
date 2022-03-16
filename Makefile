all: mat_a.cpp	
	g++ mat_a.cpp -I /usr/local/include/eigen3 -std=c++14
fn: mat_fn.cpp
    g++ mat_fn.cpp -I /usr/local/include/eigen3 -std=c++14