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
