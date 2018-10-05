DIR=writing
MANU_NEW=manuscript_revised.tex
MANU_OLD=manuscript_old.tex
MANU_CHG=manuscript_changes.tex
SI_NEW=supporting_information_revised.tex
SI_OLD=supporting_information_old.tex
SI_CHG=supporting_information_change.tex


all: manuscript SI change_manuscript change_SI clean

manuscript: $(MANU_NEW)
	latexmk -f -pdf -quiet -pdflatex='xelatex -interaction=nonstopmode' $<

SI: $(SI_NEW)
	latexmk -f -pdf -quiet -pdflatex='xelatex -interaction=nonstopmode' $<

change_manuscript: $(MANU_NEW) ref.bib
	latexdiff --exclude-textcmd="section,subsection" --config="PICTUREENV=(?:figure|section|DIFnomarkup)[*]*" --graphics-markup=0 --disable-citation-markup  $(MANU_OLD) $(MANU_NEW) > $(MANU_CHG)
	latexmk -pdf -f -quiet -pdflatex="xelatex -interaction=nonstopmode" $(MANU_CHG)

change_SI: $(SI_NEW) ref.bib
	latexdiff --exclude-textcmd="section,subsection" --config="PICTUREENV=(?:figure|section|DIFnomarkup)[*]*" --graphics-markup=0 --disable-citation-markup  $(SI_OLD) $(SI_NEW) > $(SI_CHG)
	latexmk -pdf -f -quiet -pdflatex="xelatex -interaction=nonstopmode" $(SI_CHG)

clean:
	latexmk -c

cleanall:
	latexmk -C
	rm change*

