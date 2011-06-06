SHELL = /bin/sh
UNAME = $(shell uname)

# Root variables
ROOTCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS     = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS    = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs)

################
##linux
CXX          = g++
RM           = rm -f
AR           = ar cru

##Flags
CXXFLAGS     = -O0 -g -fPIC -funroll-loops -std=c++0x -Wall

EXENAME		= imagiro
SRCEXT   	= cpp
SRCDIR  	= src
UNFOLDINGSRCDIR	= unfolding/src
INCDIR   	= include
UNFOLDINGINCDIR	= unfolding/include
OBJDIR   	= build
UNFOLDINGOBJDIR	= unfolding/build
EXEDIR  	= bin
SRCS    	:= $(shell find $(SRCDIR) -name '*.$(SRCEXT)')
UNFOLDINGSRCS 	:= $(shell find $(UNFOLDINGSRCDIR) -name '*.$(SRCEXT)')
OBJS    	:= $(patsubst $(SRCDIR)/%.$(SRCEXT),$(OBJDIR)/%.o,$(SRCS))
UNFOLDINGOBJS 	:= $(patsubst $(UNFOLDINGSRCDIR)/%.$(SRCEXT),$(UNFOLDINGOBJDIR)/%.o,$(UNFOLDINGSRCS))

GARBAGE  = $(OBJDIR)/*.o $(UNFOLDINGOBJDIR)/*.o $(EXEDIR)/$(EXENAME)

#################
##Dependencies
# Linux
ifeq "$(UNAME)" "Linux"
RANLIB       = ranlib
CXXFLAGS    += -I$(INCDIR) -I$(UNFOLDINGINCDIR) $(ROOTCFLAGS) #-I$(GSLINC)
LINKFLAGS    = -g $(shell root-config --nonew) $(shell root-config --ldflags)
endif

# OS X
ifeq "$(UNAME)" "Darwin"
RANLIB       = ranlib
CXXFLAGS    += -I$(INCDIR) -I$(UNFOLDINGINCDIR) $(ROOTCFLAGS) #-I$(GSLINC)
LINKFLAGS    =
endif

##Libraries
LIBS       += $(ROOTLIBS) -lHtml -lThread

##Targets
all : $(EXEDIR)/$(EXENAME)

$(EXEDIR)/$(EXENAME) : $(OBJS) $(UNFOLDINGOBJS)
	$(CXX) -o $@ $(OBJS) $(UNFOLDINGOBJS) $(LINKFLAGS) $(LIBS)

$(OBJDIR)/%.o : $(SRCDIR)/%.$(SRCEXT)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(UNFOLDINGOBJDIR)/%.o : $(UNFOLDINGSRCDIR)/%.$(SRCEXT)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean   :
	$(RM) $(GARBAGE)

cleanall:
	$(RM) $(GARBAGE)
