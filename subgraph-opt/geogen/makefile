ANN=./ann_1.1.2
BOOST=/usr/include/boost #../boost_1_61_0
COMP=-I${ANN}/include -L${ANN}/lib -lANN
COMP1=-I${BOOST}

geogen: geogen.c	  
	g++ -o geogen -O2 geogen.c ${COMP} ${COMP1}

debug: geogen.c
	g++ -o geogen -g geogen.c ${COMP} ${COMP1}
