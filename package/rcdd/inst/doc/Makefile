
all : design.pdf vinny.pdf clean

design.pdf : design.tex
	texi2dvi --pdf design.tex

vinny.pdf : vinny.tex
	texi2dvi --pdf vinny.tex

clean :
	rm -f *.dvi *.aux *.log vinny-fig* Rplots.ps
