CPP 		= g++
CXXFLAGS	= -g -O3 -Wall -fPIC -D_REENTRANT -Wno-deprecated -fpermissive -std=c++11

ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTLIBS     	:= $(shell root-config --libs)
ROOTGLIBS    	:= $(shell root-config --glibs)
CXXFLAGS	+= $(ROOTCFLAGS)

LIBS 		= $(ROOTGLIBS) -lMinuit $(WCSIM_BUILD_DIR)/lib/libWCSimRoot.so

INC = $(WCSIM_BUILD_DIR)/include/WCSim
SRC= $(WCSIM_BUILD_DIR)/src

CXXFLAGS += -I$(LIBS) -I$(SRC) -I$(INC) 

TARGET=   LEAFOutputAnalysisHybrid_leafclass_new #ProducePDF_new ProducePDF_DIR   FitVertexLE AnalyzeWSHierarchy LEAFOutputAnalysisHybrid_master LEAFOutputAnalysisHybrid_leafclass

all: $(TARGET)
LEAFOutputAnalysisHybrid_leafclass: LEAFOutputAnalysisHybrid_leafclass.o
LEAFOutputAnalysisHybrid_leafclass_new: LEAFOutputAnalysisHybrid_leafclass_new.o
LEAFOutputAnalysisHybrid_master: LEAFOutputAnalysisHybrid_master.o
ProducePDF: ProducePDF.o
ProducePDF_new: ProducePDF_new.o
ProducePDF_DIR: ProducePDF_DIR.o
AnalyzeWSHierarchy: AnalyzeWSHierarchy.o
FitVertexLE: FitVertexLE.o

%: %.o
	@echo "Now make $@"
	@$(CPP) -o $@ $< $(CXXFLAGS) $(LIBS) 
	@echo "..Compile done! "

%.o: %.c
	@echo "$<"
	@echo "Start Compiling $<"
	@$(CPP) $(CXXFLAGS) -c $<
	@echo ".. Compiling Object Files $<   --> done"
	@echo ""

%.o: %.cpp
	@echo "$<"
	@echo "Start Compiling $<"
	@$(CPP) $(CXXFLAGS) -c $<
	@echo ".. Compiling Object Files $<   --> done"
	@echo ""

%.o: %.C %.h
	@echo '<< compiling' $< '>>'
	@$(CXX) $(CXXFLAGS) $(INCFLAGS) -c -o $@ $<

clean: 
	@echo "Now Clean Up"
	rm -f $(TARGET) *~ *.o *.o~ core
