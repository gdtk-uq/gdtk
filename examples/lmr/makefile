test:
	PYTHONDONTWRITEBYTECODE=1 pytest -v

examples-listing.adoc: examples-listing-control.adoc
	asciidoctor-reducer -o $@ examples-listing-control.adoc

examples-listing.html: examples-listing.adoc
	asciidoctor examples-listing.adoc

clean:
	- rm examples-listing.adoc
	- rm examples-listing.html

