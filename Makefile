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

.PHONY: test
test:
	cd test && \
		python AXIsolver.py
