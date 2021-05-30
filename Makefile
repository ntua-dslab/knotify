.PHONY: venv grammars clean clean-yaep clean-venv deps

all: grammars venv

#####################################################
# Dependencies

deps:
	sudo apt-get install virtualenv git gcc make g++ build-essential bison libgsl23

#####################################################
# YAEP

YAEP_DIR = .yaep
CFLAGS += -I$(YAEP_DIR)/src
LIBS += $(YAEP_DIR)/src/libyaep.a

grammars: libpseudoknot.so libhairpin.so

lib%.so: yaep_parsers/%.c $(YAEP_DIR)/src/libyaep.a
	$(CC) $< $(CFLAGS) $(LIBS) -fPIC -shared -o $@

$(YAEP_DIR)/src/libyaep.a:
	git clone https://github.com/ntua-dslab/yaep --depth 1 $(YAEP_DIR)
	cd $(YAEP_DIR) && ./configure CFLAGS=-fPIC && make -j

#####################################################
# Python

VENV_DIR = .venv

venv: $(VENV_DIR)/bin/rna_analysis

$(VENV_DIR)/bin/rna_analysis: setup.py setup.cfg
	virtualenv --python=python3 $(VENV_DIR)
	$(VENV_DIR)/bin/pip install -e . -r wheel-requirements.txt -r dev-requirements.txt

#####################################################
# Clean

clean: clean-yaep clean-venv

clean-yaep:
	rm -rf $(YAEP_DIR) **.so

clean-venv:
	rm -rf $(VENV_DIR)
