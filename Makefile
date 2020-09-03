paths: paths.cpp
	g++ -std=c++14 $^ -lSDL2 -lGL -o $@
