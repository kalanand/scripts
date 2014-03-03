# Please see SoftRelTools/HOWTO-GNUmakefile for documentation
# $Id: GNUmakefile,v 1.4 2006/05/19 17:34:18 fwinkl Exp $
#################################################################
#++ library products				[build it with 'lib']

LIBREMOVEFILES :=
LIBTMPLFILES :=
LIBDDLORDERED :=

#################################################################
#++ extra binary products	[not in production, build it with extrabin]

EXTRABINS := genShapeSyst genShapeSyst-manual BdkProfilePdf

$(addprefix $(bindir),$(EXTRABINS)): $(bindir)% : %.o

#################################################################
#++ binary products				[build it with 'bin']

BINS :=
BINCCFILES := $(BINS:=.cc) $(EXTRABINS:=.cc)

#++ Binary rules		 [in production, build it with 'bin']

$(addprefix $(bindir),$(BINS)): $(bindir)% : %.o

#++ shell script products.. 			[build it with 'bin']
BINSCRIPTS := 

#################################################################
#++ regression test scripts			[build it with 'test']

# $(testdir)mytest.T : mytest.tcl mytesttemp2

#################################################################
#++ include standard makefile from SoftRelTools.
include SoftRelTools/standard.mk
