TXT = \
    cnv-calling.txt \
    consensus-sequence.txt \
    install.txt \
    index.txt \

HTML = $(TXT:.txt=.html)
XML  = $(TXT:.txt=.xml)

all: $(HTML)

%.html: %.txt header.txt footer.txt
	asciidoctor -b html5 -a linkcss -a stylesheet=index.css -a no-header-footer $<

clean:
	@rm -f $(HTML) $(XML)

