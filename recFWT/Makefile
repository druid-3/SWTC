################################################################################
################################################################################
##  developed by DRUID's labs.                                                ##
##  druid3@i.ua    druid-iii@yandex.ru                                        ##
##  http://www.druid3.cc.ua                                                   ##
##                                                                            ##
##                 00  000   0  0  0  000   0   000                           ##
##                000  0  0  0  0  0  0  0     00        *******              ##
##               0100  000   0  0  0  0  0       00     * ****  *             ##
##              01100  0  0   000  0  000      000     *  *   *  *            ##
##             011 00                                  *  ****   *            ##
##            011  00     11   1      111               * *   * *             ##
##           00000000    1  1  111   11                  *******              ##
##          00000000011  1111  1  1    11                                     ##
##          11111111111  1  1  111   111                                      ##
##                                                                            ##
################################################################################
################################################################################
# Clean automated Makefile 4 C&C++ (prg) 
# distributed under GNU GPL v.3
WHERE_IS_GCC  = /usr/bin/

#INCLUDE=-Isrc/include -I../src/include

CC            = $(WHERE_IS_GCC)gcc
CPP           = $(WHERE_IS_GCC)g++
STRIP         = $(WHERE_IS_GCC)strip

CC_FLAGS      = -O2                                                            \
                -W                                                             \
                -Wall                                                          \
                -Winline                                                       \
                -Wmissing-noreturn                                             \
                -Wmissing-prototypes                                           \
                -Wnested-externs                                               \
                -Wconversion                                                   \
                -Wfloat-equal                                                  \
                -Wshadow                                                       \
                -Wsequence-point                                               \
                -pipe 
# добавить, если разрабатывается модуль ядра 
# <add, if you are developing a kernel module>
#KERNEL        = -DMODULE -D__KERNEL__ -DLINUX
#INCLUDE=-I/usr/include/linux/

LD_FLAGS      = -L./LIB/ $(ALIB) $(DLIB)
DLIB          =
ALIB          =


PRG_NAME      = universalRecursiveFWT_TEST

PRJDIRS       = ./OBJ                                                          \
                ./ASM                                                          \
                ./LIB

OBJFILES     := $(patsubst %.c,   ./OBJ/%.o,       $(wildcard *.c))
OBJFILES     += $(patsubst %.cpp, ./OBJ/%.o,       $(wildcard *.cpp))
OBJFILES     += $(patsubst %.s,   ./OBJ/%.o,       $(wildcard *.s))

2ASM         := $(patsubst %.c,   ./ASM/%.S,       $(wildcard *.c))
2ASM         += $(patsubst %.cpp, ./ASM/%.S,       $(wildcard *.cpp))


#  include <bla-bla-bla> предупреждение и останов если нет Makefile 
# -include <bla-bla-bla> нет предупреждения и останова если нет Makefile 
-include ./LIB/Makefile

all: prj_dirs myprog show_src run_script

#PRG----------------------------------------------------------------------------
myprog: $(OBJFILES) 
	$(CC) -o $(PRG_NAME) $(OBJFILES) $(LD_FLAGS)

./OBJ/%.o: %.c
	$(CC) $(CC_FLAGS) -c -o $@ $<
./OBJ/%.o: %.cpp
	$(CPP) $(CC_FLAGS) -c -o $@ $<
./OBJ/%.o: %.s  
	$(CC) -o $@ -c $<
#-------------------------------------------------------------------------------

run: all
	./run

src2asm: $(2ASM)
./ASM/%.S: %.[cpp,c] 
	$(CC) -o $@ -S -fverbose-asm $<

prj_dirs:
	@$(foreach DIR, $(PRJDIRS),                                            \
	if ! test -d $(DIR);                                                   \
		then                                                           \
			mkdir $(DIR);                                          \
			chmod 775 $(DIR);                                      \
		else                                                           \
			echo ' папка $(DIR) существует ';                      \
	fi;                                                                    \
	)

run_script:
	@if ! test -f ./run;                                                   \
	then                                                                   \
		touch ./run;                                                   \
		chmod 775 ./run;                                               \
		echo '#!/bin/bash'                                   >> ./run; \
		echo 'echo "starting ..."'                           >> ./run; \
		echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./LIB' >> ./run; \
		echo 'exec ./$(PRG_NAME) $$@'                        >> ./run; \
	else                                                                   \
		echo 'файл run существует';                                    \
	fi

#-------------------------------------------------------------------------------
# * !!! внимание !!! 
# к проекту присоединяются ВСЕ исходные тексты из данной директории
# <ALL the source code from this directory join to the project.>
# * !!! <attention> !!!
show_src:
	@echo "-------------------------------------------------------------------------------"
	@echo "src in project:"
	@echo
	@for i in `ls *.[c,h,hpp,cpp,s]`;                                                     \
		do                                                                            \
			echo $${i};                                                           \
		done
	@echo "-------------------------------------------------------------------------------"
#-------------------------------------------------------------------------------
count:
	@wc -l *.[c,h,cpp,hpp,s]

archive: clean
	@tar -zcvf ../$(PRG_NAME)__`date +%y.%m.%d`.tgz ./ --recursion -P

clean:
	@rm -f $(PRG_NAME) run
	@$(foreach DIR, ./ $(PRJDIRS), rm -f                                   \
	                                     $(DIR)/*.o                        \
	                                     $(DIR)/*.d                        \
	                                     $(DIR)/*.map                      \
	                                     $(DIR)/*.hex                      \
	                                     $(DIR)/*.listing                  \
	                                     $(DIR)/*.bak                      \
	                                     $(DIR)/*.out                      \
	                                     $(DIR)/*.swp                      \
	                                     $(DIR)/*~                         \
	                                     $(DIR)/core                       \
	                                     $(DIR)/*.tu                       \
	  )

striping:
	@$(STRIP) ./$(PRG_NAME)

doc++:
	@doc++ --all --html --tables-border --private --dir ./DOC *.[c,h,hpp,cpp,s]

.PHONY: all
.PHONY: myprog
.PHONY: run
.PHONY: src2asm
.PHONY: prj_dirs
.PHONY: run_script
.PHONY: show_src
.PHONY: count
.PHONY: archive
.PHONY: clean
.PHONY: striping
.PHONY: doc++


