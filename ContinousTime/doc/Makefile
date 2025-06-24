##
##  Simple makefile to create a paper
##

TARGET= ContinousTimeAlgorithm

paper:
	rm *.dvi *.out *.toc *.tmp *.log *.aux *.bbl *.blg;\
	pdflatex $(TARGET).tex;
	bibtex $(TARGET);\
	pdflatex $(TARGET).tex;\
	pdflatex $(TARGET).tex;\
#	dvips -t letter -Ppdf -z -G0  $(TARGET).dvi -xo
#	ps2pdf  $(TARGET).ps

FIG = fig

BIB = bib


ALL_FILES = Makefile $(TARGET).tex $(TARGET).pdf fig/
# bib/ fig/ nb/

tar: $(ALL_FILES) 
	tar -cvjf $(TARGET).tar $(ALL_FILES) 

clean:
	rm *.dvi *.out *.toc *.tmp *.log *.aux *.bbl *.blg *.*~

figures.tar: $(FIGURES) 
	tar cvf figures.tar $(FIGURES) 


