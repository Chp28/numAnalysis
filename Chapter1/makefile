run: B.cpp C.cpp D.cpp E.cpp F.cpp
	g++ B.cpp -o B
	./B
	g++ C.cpp -o C
	./C
	g++ D.cpp -o D
	./D
	g++ E.cpp -o E
	./E
	g++ F.cpp -o F
	./F
	
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
