v.png: plotter.py output.txt
	python plotter.py
	rm -f output.txt placas

output.txt: placas
	./placas > output.txt

placas: placas.c
	gcc -o placas placas.c
