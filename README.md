# APLCON C++11 wrapper

APLCON++ is a C++11 wrapper class around the Constrained Least Squares
fitter APLCON by Volker Blobel. Please <a href="http://www.desy.de/~blobel/wwwcondl.html">see
his DESY homepage how APLCON works</a>. Also, if you use this wrapper
in a publication, don't forget to cite Volker Blobel's original work:

```latex
\bibitem{bib:blolo}
Volker Blobel and Erich Lohrmann,
\emph{Statistische Methoden der Datenanalyse},
Teubner Studienb\"{u}cher, Teubner (1998);
e-book \url{http://www.desy.de/~blobel/eBuch.pdf}
```

Please report bugs and improvements!

## Installation and Usage

To build this project, you need
  * GNU compiler version >4.9.2 with g++ and gfortran
  * cmake version >3.0

Then do a standard out-of-source build by executing the following
inside the projects directory:

    mkdir build
    cd build
    cmake ..
    make

You'll find a library called `libaplcon++.a` inside the above created
`build` directory, the corresponding header to include is
`src/APLCON.hpp`.

## Tests and benchmarking

You can build benchmarking binaries by calling `cmake
-DEnableBenchmark=On ..` and running `./bench/bench_APLCON`. Tests are
build and run by `make build_and_test`.
