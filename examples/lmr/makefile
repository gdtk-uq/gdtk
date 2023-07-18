test:
	PYTHONDONTWRITEBYTECODE=1 pytest -v

examples-listing.adoc: elc.adoc
	asciidoctor-reducer -o $@ $<

examples-listing.html: examples-listing.adoc
	asciidoctor $<

clean:
	- rm examples-listing.adoc
	- rm examples-listing.html

