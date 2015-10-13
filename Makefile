JULIA=julia

FILES = $(shell ls plots)

all:
	for file in $(FILES); do \
		$(JULIA) plots/$$file; done
