JULIA=julia

FILES = $(shell find plots tables)

all:
	for file in $(FILES); do \
		$(JULIA) $$file; done
