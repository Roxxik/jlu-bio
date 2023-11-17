.PHONY: all clean

DIAMOND = /vol/software/bin/diamond-2.0.0

# overwrite by calling `make all TARGET=foo`
TARGET=Bacillus

# Get all the .faa files from the proteomes directory except the target
OTHER_FAA_FILES = $(filter-out proteomes/$(TARGET).faa, $(wildcard proteomes/*.faa))

# Extract the base names of the other .faa files
OTHER_BASE_NAMES = $(notdir $(basename $(OTHER_FAA_FILES)))

# create an empty list
TSV_FILES :=

define BLAST_RULE_TEMPLATE
TSV_FILES += out/$(TARGET)_v_$(1).tsv out/$(1)_v_$(TARGET).tsv

out/$(1)_v_$(TARGET).tsv: out/$(TARGET).dmnd proteomes/$(1).faa
	$$(DIAMOND) blastp -d out/$(TARGET).dmnd -q proteomes/$(1).faa -o out/$(1)_v_$(TARGET).tsv

out/$(TARGET)_v_$(1).tsv: out/$(1).dmnd proteomes/$(TARGET).faa
	$$(DIAMOND) blastp -d out/$(1).dmnd -q proteomes/$(TARGET).faa -o out/$(TARGET)_v_$(1).tsv
endef

# Dynamically generate rules to compare target with other .faa files (both ways) and populate TSV_FILES
$(foreach subject, $(OTHER_BASE_NAMES), \
    $(eval $(call BLAST_RULE_TEMPLATE,$(subject))) \
)

# Rules to create diamond databases from .faa files in proteomes
out/%.dmnd: proteomes/%.faa
	$(DIAMOND) makedb --in $< -d $@

all: out/$(TARGET)_reciprocal.tsv

out/$(TARGET)_reciprocal.tsv: script2.py $(TSV_FILES)
	python script2.py proteomes/$(TARGET).faa -o $@

clean:
	rm -f out/*.tsv
	rm -f out/*.dmnd