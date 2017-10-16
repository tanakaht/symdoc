all: sqrt.html

clean:
	rm sqrt.md sqrt.html
	rm -rf __pycache__

sqrt.html: sqrt.py
	python sqrt.py --symdoc
	pandoc -t html -o sqrt.html --standalone --mathjax sqrt.md