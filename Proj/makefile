run:
	g++ -I ./eigen main.cpp -o main
	./main
report:
	latex design.tex
	dvipdf design.dvi
	rm design.log
	rm design.aux
	rm design.dvi
	latex report.tex
	dvipdf report.dvi
	rm report.log
	rm report.aux
	rm report.dvi
