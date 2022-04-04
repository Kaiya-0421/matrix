all: mat_a.cpp	
	g++ mat_a.cpp -I /usr/local/include/eigen3 -std=c++14
fn: mat_fn.cpp
	g++ mat_fn.cpp -I /usr/local/include/eigen3 -std=c++14
mat: mat.cpp
	g++ mat.cpp -I /usr/local/include/eigen3 -std=c++14
bl: mat_ch_bl.cpp
	g++ mat_ch_bl.cpp -I /usr/local/include/eigen3 -std=c++14 