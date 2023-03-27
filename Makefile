DIR_TOPLEVEL := $(shell git rev-parse --show-toplevel)

all: build/Makefile
	$(MAKE) -C build all

build/Makefile:
	cmake -B build -DCMAKE_BUILD_TYPE=RELEASE

clean:
	rm -rf build/
