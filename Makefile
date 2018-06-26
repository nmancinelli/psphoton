all: bin/psphoton bin/post

bin/psphoton: src/psphoton.f90 Makefile
	gfortran -O src/psphoton.f90 -o bin/psphoton

bin/post: src/post.f90 Makefile
	gfortran src/post.f90 -o bin/post

.PHONY: clean

clean:
	rm bin/*
