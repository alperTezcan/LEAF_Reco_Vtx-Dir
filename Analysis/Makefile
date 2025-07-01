CPP 		= g++
CXXFLAGS	= -g -O3 -Wall -fPIC -D_REENTRANT -Wno-deprecated -fpermissive -std=c++11

ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTLIBS     	:= $(shell root-config --libs)
ROOTGLIBS    	:= $(shell root-config --glibs)
CXXFLAGS	+= $(ROOTCFLAGS)

LIBS 		= $(ROOTLIBS) $(ROOTGLIBS) $(WCSIM_BUILD_DIR)/lib/libWCSimRoot.so
#LIBS 		= $(ROOTLIBS) $(ROOTGLIBS) -Wl,-rpath,$(WCSIM_BUILD_DIR)/libWCSimRoot.so

INC = $(WCSIM_BUILD_DIR)/include/WCSim
SRC= $(WCSIM_BUILD_DIR)/src

CXXFLAGS += -I$(SRC) -I$(INC) 
#CXXFLAGS += -I$(LIBS) -I$(SRC) -I$(INC) 

TARGET= CompareProfiles #dtowall_info # AnalyzeWSHierarchy_new  #AnalyzeWSHierarchy AnalyzeFQFull

all: $(TARGET)
AnalyzeWSHierarchy: AnalyzeWSHierarchy.o
AnalyzeWSHierarchy_new: AnalyzeWSHierarchy_new.o
AnalyzeFQFull: AnalyzeFQFull.o
CompareProfiles: CompareProfiles.o
dtowall_info: dtowall_info.o

%: %.o
	@echo "Now make $@"
	@echo $(LIBS)
	@$(CPP) -o $@ $< $(LIBS)
#@$(CPP) -o $@ $< $(CXXFLAGS) $(LIBS)
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

clean: 
	@echo "Now Clean Up"
	rm -f $(TARGET) *~ *.o *.o~ core
