SRCS = $(wildcard src/*.f)
OBJS = $(SRCS:.f=.o)
OBJS += src/clear.o
FOPTS = --std=legacy -ffpe-trap=invalid,zero
FC = gfortran


all: $(OBJS)
	$(FC) $(FOPTS) -o rrbaxi $^

%.o:	%.f
	$(FC) $(FOPTS) -c -o $@ $<

src/clear.o:	src/clear.c
	gcc -c -o $@ $<

clean:
	rm -f src/.o rrbaxi

.PHONY: venv
venv:
	echo 'layout python3' > .envrc && \
		direnv allow

.PHONY: init
init:
	pip install -U pip
	pip install pip-tools

.PHONY: reqs
reqs:
	pip-compile
	pip install -r requirements.txt
	jupyter contrib nbextensions install
	cp ~/_sys/tikzmagic.py .direnv/python-3.10.8/lib/python3.10/site-packages

.PHONY: nb
nb:
	cd  book && \
		jupyter notebook

.PHONY: book
book:
	jb build book
	cp -R book/_build/html/* docs


.PHONY: test
test:
	cd test && \
		python AXIsolver.py
