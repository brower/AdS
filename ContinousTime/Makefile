.SUFFIXES:
.SUFFIXES: .o .cpp
#============================================================



TARGET1 = ConTime
C_OBJS1 = ConTime.o



C_SOURCES = ConTime.cpp
MY_INCLUDES = 

ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES) plotMult.py ./doc ./ref

CCX = g++ 
CXXFLAGS = -g -Wall
#-std=c99

#============================================================
all: $(TARGET1) 

.o:.cpp	$(MY_INCLUDES)
	$(CCX)  -c  $(CXXFLAGS) $<  

$(TARGET1) :   $(C_OBJS1)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@

$(TARGET2) :   $(C_OBJS2)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@



# Implicit rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites
#============================================================


clean:
	rm -f $(TARGET1) $(TARGET2) $(TARGET3)  $(C_OBJS1) $(C_OBJS2) $(C_OBJSTEST1) $(C_OBJ3) *~

tar: $(ALL_SOURCES)
	tar cvf $(TARGET1).tar $(ALL_SOURCES)



