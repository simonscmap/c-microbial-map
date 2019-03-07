SINGULARITY = /usr/local/bin/singularity
IMG = c-microbial-map-0.0.1.img
DEF = image.def

run:
	sudo $(SINGULARITY) run $(IMG)

clean:
	rm -f $(IMG)

img: clean
	sudo $(SINGULARITY) build $(IMG) $(DEF)

img_old: clean
	sudo $(SINGULARITY) create --size 2048 $(IMG)
	sudo $(SINGULARITY) bootstrap $(IMG) $(DEF)

shell:
	sudo $(SINGULARITY) shell --writable -B $(shell pwd):/home $(IMG)
