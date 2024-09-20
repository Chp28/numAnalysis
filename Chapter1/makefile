B: B.cpp
	g++ B.cpp -o B
	g++ C.cpp -o C
	g++ D.cpp -o D
	g++ E.cpp -o E
	g++ F.cpp -o F
	
report: answer.tex
	latex answer.tex
	dvipdf answer.dvi
	rm answer.log
	rm answer.aux
	rm answer.dvi
	
clean:
	rm B
	rm C
	rm D
	rm E
	rm F
