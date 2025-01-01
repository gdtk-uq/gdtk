# Look for pytest under different names
PYTEST := $(shell which pytest)
PYTEST := $(or $(PYTEST), pytest-3)

test:
	PYTHONDONTWRITEBYTECODE=1 $(PYTEST) -v

examples-listing.adoc: elc.adoc
	asciidoctor-reducer $< | awk -f insert-imagesdir.awk > $@

examples-listing.html: examples-listing.adoc
	asciidoctor $<

clean:
	- rm examples-listing.adoc
	- rm examples-listing.html

