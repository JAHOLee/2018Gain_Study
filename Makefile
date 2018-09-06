CXXFLAGS=-g -m64 -O2 -Wall -std=c++0x
ROOTFLAGS=$(shell root-config --libs --cflags --glibs)

makePlots: main.o makePlots.o fitter.o compare.o single_mod.o
	g++ $^ -o makePlots $(CXXFLAGS) $(ROOTFLAGS)
	mv *.o obj_file

main.o: main.cc makePlots.h fitter.h compare.h
	g++ -c $(CXXFLAGS) $(ROOTFLAGS) $< -o $@

makePlots.o:makePlots.cc makePlots.h fitter.h
	g++ -c $(CXXFLAGS) $(ROOTFLAGS) $< -o $@

fitter.o:fitter.cc fitter.h
	g++ -c $(CXXFLAGS) $(ROOTFLAGS) $< -o $@

compare.o:compare.cc compare.h
	g++ -c $(CXXFLAGS) $(ROOTFLAGS) $< -o $@

single_mod.o:single_module.cc single_module.h
	g++ -c $(CXXFLAGS) $(ROOTFLAGS) $< -o $@

clean:
	rm -f makePlots obj_file/*.o ./*~
