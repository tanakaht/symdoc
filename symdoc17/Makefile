all: kk3.html

clean:
	rm kk3.md kk3.html
	rm -rf __pycache__

kk3.html: kk3.py
	python kk3.py --symdoc
	pandoc -t html -o kk3.html --standalone --mathjax kk3.md