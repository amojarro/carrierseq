run: clean build

clean:
	rm -f cseq

build: clean
	sudo singularity build cseq Singularity
