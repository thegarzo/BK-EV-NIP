FFTW_INC=/usr/local/include
FFTW_DIR=/usr/local/lib

CXX           = g++
CXXFLAGS      = -std=c++17 -Wall -fPIC -O3 -I/opt/homebrew/include/eigen3
LD            = g++
LDFLAGS       = -O3

LIBS          =  -I/usr/local/include -L/usr/local/lib -I${FFTW_INC} -L${FFTW_DIR} -lgsl -lgslcblas -lcuba -lm `lhapdf-config --cflags --ldflags` -lfftw3

BUILD = $(PWD)/build
MCDIP = $(PWD)

vpath %.cpp src
objdir     = obj
tabdir     = tabs

SRC        = main.cpp dipole.cpp
OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC))

TARGET	   = compute
#------------------------------------------------------------------------------
$(TARGET): $(OBJS) $(TABS)
		$(LD)  $(LDFLAGS) $^ -o $@ $(LIBS)
		@echo "$@ done"


$(OBJS): | $(objdir)

$(objdir):
	@mkdir -p $(objdir)

$(TABS): | $(tabdir)

$(tabdir):
	@mkdir -p $(tabdir)

obj/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

builder:
		if [ ! -d $(BUILD) ]; then make buildfldr; fi
		cd $(BUILD)



clean:
		@rm -f $(OBJS) $(TARGET)

again:
		make clean
		make

renew:
	clear
	make

