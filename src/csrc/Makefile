# xtea-scan build. Links htslib via pkg-config — NO hardcoded paths (works for local, conda, container).
#
# Local build: point pkg-config at the htslib 1.22 install, e.g.
#   PKG_CONFIG_PATH=/home/cos689/software/lib/pkgconfig make
# In conda/container the packaged htslib's .pc is already on PKG_CONFIG_PATH.

CXX     ?= g++
CXXFLAGS ?= -O2 -Wall -Wextra -std=c++17 -pthread
PKG     ?= pkg-config
comma   := ,

HTS_CFLAGS := $(shell $(PKG) --cflags htslib)
HTS_LIBS   := $(shell $(PKG) --libs htslib)
HTS_LIBDIR := $(shell $(PKG) --variable=libdir htslib)
# Embed an rpath so the binary finds libhts at runtime without LD_LIBRARY_PATH.
RPATH      := $(if $(HTS_LIBDIR),-Wl$(comma)-rpath$(comma)$(HTS_LIBDIR),)

BIN := xtea_scan
SRC := xtea_scan.cpp

.PHONY: all clean check-hts
all: check-hts $(BIN)

check-hts:
	@$(PKG) --exists htslib || { echo "ERROR: htslib not found by pkg-config. Set PKG_CONFIG_PATH."; exit 1; }
	@echo "htslib $$($(PKG) --modversion htslib) via $(HTS_LIBDIR)"

$(BIN): $(SRC)
	$(CXX) $(CXXFLAGS) $(HTS_CFLAGS) -o $@ $< $(HTS_LIBS) $(RPATH)

clean:
	rm -f $(BIN)
