# set "EISPACK= /home/iit/Documents/General6R/EISPACK"
# gfortran -c *
# cd IK_CODE
# gfortran -c main.c matrix.c read.f read_ip.c eigen1.c improve.c reduce.c degen.c setup.c setup2.c check_degen.c gauss.c
# gfortran -static -o General6R main.o matrix.o read.o read_ip.o eigen1.o improve.o reduce.o degen.o setup.o setup2.o check_degen.o gauss.o %EISPACK%/dgeco.o %EISPACK%/dgedi.o %EISPACK%/rg.o %EISPACK%/rgg.o %EISPACK%/hybrj.o %EISPACK%/dsvdc.o %EISPACK%/dgesl.o

# C:\Program Files (x86)\GnuWin32\bin
# C:\MinGW\bin

#DEFS = -DIK_OPT
CC = gcc
CFLAGS = -O ${DEFS}
FFLAGS = 
LDFLAGS =  
COMPILER = gfortran
# this should be set to the corresponding directory
EISPACK = "$/home/iit/6R_G/EISPACK"
IK_CODE = "$/home/iit/6R_G/IK_CODE"
CFILES = *.c

IK_object = main.o matrix.o read.o read_ip.o eigen1.o improve.o reduce.o \
	degen.o setup.o setup2.o check_degen.o gauss.o

	
EISCODE = $(EISPACK)/dgeco.o $(EISPACK)/dgedi.o $(EISPACK)/rg.o $(EISPACK)/rgg.o $(EISPACK)/hybrj.o $(EISPACK)/dsvdc.o $(EISPACK)/dgesl.o

#.f.o : ;cd $(EISPACK); $(COMPILER) -c $(FFLAGS) $<; cd $(IK_CODE)
.f.o : ;$(COMPILER) -c $(FFLAGS) $<

.c.o : ; $(CC) -c $(CFLAGS) $<

LIB= eiscode
OJBS =test

scratch:
	make eiscode
	make test

test : $(IK_object)
	$(COMPILER) $(CFLAGS) -o test1 $(IK_object) $(EISCODE)


eiscode: 
	cd $(EISPACK); $(COMPILER) -c *; cd $(IK_CODE)

clean: 
	rm *.o


depend :
	makedepend -- $(CFLAGS) -- $(CFILES)

# DO NOT DELETE THIS LINE -- make depend depends on it.
