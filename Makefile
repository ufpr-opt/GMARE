JULIA=julia

FILES = $(shell ls plots tables)

all:
	for file in $(FILES); do \
		$(JULIA) plots/$$file; done
