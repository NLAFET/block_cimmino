sed -n '
1h
1!H
$ {
  g
  s/\\begin{fulllineitems}\n*\\end{fulllineitems}//g
  p
}' build/latex/ABCDSolver.tex
