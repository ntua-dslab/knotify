.PHONY: venv grammars clean clean-yaep clean-venv deps pairaligns

all: libraries venv

libraries: grammars energies pairaligns

#####################################################
# Dependencies

deps:
	sudo apt-get update
	sudo apt-get install -y virtualenv git gcc make g++ build-essential bison libgsl23

#####################################################
# YAEP

YAEP_DIR = .yaep
CFLAGS += -I$(YAEP_DIR)/src
LIBS += $(YAEP_DIR)/src/libyaep.a

grammars: libpseudoknot.so libhairpin.so libbruteforce.so

lib%.so: parsers/%.c $(YAEP_DIR)/src/libyaep.a
	$(CC) $< $(CFLAGS) $(LIBS) -fPIC -shared -o $@

$(YAEP_DIR)/src/libyaep.a:
	git clone https://github.com/ntua-dslab/yaep --depth 1 $(YAEP_DIR)
	cd $(YAEP_DIR) && ./configure CFLAGS=-fPIC && make

#####################################################
# PK energy

PKENERGY_DIR = pkenergy

energies: libpkenergy.so

libpkenergy.so:
	cd $(PKENERGY_DIR)/hotknots/LE && make -j
	cp $(PKENERGY_DIR)/hotknots/LE/libpkenergy.so .

#####################################################
# PairAligns

pairaligns: libskipfinalau.so libcpairalign.so

libskipfinalau.so: pairalign/skipfinalau.c
	$(CC) $<  -fPIC -shared -o $@

libcpairalign.so: pairalign/cpairalign.c
	$(CC) $<  -fPIC -shared -o $@

#####################################################
# Python

VENV_DIR = .venv

venv: $(VENV_DIR)/bin/rna_analysis

$(VENV_DIR)/bin/rna_analysis: setup.py setup.cfg
	virtualenv --python=python3 $(VENV_DIR)
	$(VENV_DIR)/bin/pip install -e . -r wheel-requirements.txt -r dev-requirements.txt

#####################################################
# IPKnot

IPKNOT_DIR = .ipknot

ipknot: $(IPKNOT_DIR)/ipknot

ipknot-deps:
	sudo apt-get update
	wget https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_20_04/viennarna-dev_2.5.0-1_amd64.deb -O viennarna.deb
	sudo apt-get install -y ./viennarna.deb cmake build-essential pkg-config libltdl-dev libglpk-dev liblz4-dev

$(IPKNOT_DIR)/ipknot: ipknot-deps
	git clone https://github.com/satoken/ipknot --depth 1 $(IPKNOT_DIR)
	cd $(IPKNOT_DIR) && cmake . && make -j

#####################################################
# Knotty

KNOTTY_DIR = .knotty

knotty: $(KNOTTY_DIR)/knotty

knotty-deps:
	sudo apt-get install -y cmake

$(KNOTTY_DIR)/knotty: knotty-deps
	git clone https://github.com/HosnaJabbari/Knotty --depth 1 $(KNOTTY_DIR)
	cd $(KNOTTY_DIR) && cmake . && make -j

#####################################################
# HotKnots

HOTKNOTS_DIR = .hotknots

hotknots: $(HOTKNOTS_DIR)/HotKnots_v2.0/bin/HotKnots

hotknots-deps:
	mkdir -p $(HOTKNOTS_DIR)
	sudo apt-get install -y libx11-dev

$(HOTKNOTS_DIR)/HotKnots_v2.0.tar.gz:
	mkdir -p $(HOTKNOTS_DIR)
	cd $(HOTKNOTS_DIR) && wget http://www.cs.ubc.ca/labs/beta/Software/HotKnots/HotKnots_v2.0.tar.gz && tar xvzf HotKnots_v2.0.tar.gz

$(HOTKNOTS_DIR)/HotKnots_v2.0/bin/HotKnots: hotknots-deps $(HOTKNOTS_DIR)/HotKnots_v2.0.tar.gz
	cd $(HOTKNOTS_DIR)/HotKnots_v2.0/graphics && make clean
	cd $(HOTKNOTS_DIR)/HotKnots_v2.0 && make clean
	cd $(HOTKNOTS_DIR)/HotKnots_v2.0/graphics && make -j
	cd $(HOTKNOTS_DIR)/HotKnots_v2.0 && make -j

#####################################################
# Iterative-HFold

IHFOLD_DIR = .iterative-hfold

ihfold: $(IHFOLD_DIR)/HFold_iterative

ihfold-deps:
	sudo apt-get install -y cmake

$(IHFOLD_DIR)/HFold_iterative: ihfold-deps
	git clone https://github.com/HosnaJabbari/Iterative-HFold --depth 1 $(IHFOLD_DIR)
	cd $(IHFOLD_DIR) && cmake . && make -j

#####################################################
# Clean

clean: clean-yaep clean-venv clean-libs clean-pkenergy clean-ipknot clean-knotty clean-hotknots clean-ihfold

clean-libs:
	rm -rf **.so

clean-pkenergy:
	cd $(PKENERGY_DIR)/hotknots/LE && make clean

clean-yaep:
	rm -rf $(YAEP_DIR)

clean-venv:
	rm -rf $(VENV_DIR) .pytest_cache **.egg-info .eggs

clean-ipknot:
	rm -rf $(IPKNOT_DIR)

clean-knotty:
	rm -rf $(KNOTTY_DIR)

clean-hotknots:
	rm -rf $(HOTKNOTS_DIR)

clean-ihfold:
	rm -rf $(IHFOLD_DIR)
