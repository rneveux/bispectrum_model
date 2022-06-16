all:
#	python setup.py build_ext -if --inplace && python2 setup.py build_ext -if --inplace
#	python setup.py install --user && python2 setup.py install --user
	python setup.py install --user 
clean:
	-rm -r build hitomi.cpp *.so *.pyc
