v.png: plotter.py output.txt
	python plotter.py
	rm -f output.txt placas

output.txt: Entrega
	./Entrega > output.txt

placas: Entrega.c
	mpicc -o Entrega Entrega.c
