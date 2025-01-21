This is a Rust port of the program [allcurve.pas](https://www.pdmi.ras.ru/~arnsem/dataprog/allcurve.pas). See the original copyright notice in the main source file.

The program enumerates closed curves in a plane with a given number of self-intersections. It uses the arc representation of curves, outlined in: Fedor Duzhin and Biaoshuai Tao, [On computational complexity of plane curve invariants](https://hdl.handle.net/10356/107188), Online Journal of Analytic Combinatorics, vol. 9, no. 4 (2014).

Also, this program has extra functionality.
- The oriented curves in the oriented plane are broken down by index (turning number), allowing to extend the OEIS array [A008985](https://oeis.org/A008985).
- It is checked whether the Gauss diagram of a curve has intersections, allowing to enumerate "tree-like curves" (sequence [A118814](https://oeis.org/A118814)).

Usage: there is one optional command line argument, the number of intersections.

Related programs:
- [The Java port of allcurves.pas in jOEIS](https://github.com/archmageirvine/joeis/blob/master/src/irvine/oeis/a008/A008980.java);
- J. Bétréma's [program](https://github.com/j2b2/TaitCurves) enumerating irreducible indecomposable curves;
- [plantri](https://users.cecs.anu.edu.au/~bdm/plantri/) for spherical curves.
