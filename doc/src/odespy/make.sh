#!/bin/sh

name=wrap_odespy

doconce format pdflatex $name
ptex2tex -DMINTED $name
pdflatex -shell-escape $name
exit

doconce format sphinx $name
theme=pyramid
#theme=fenics_minimal
doconce sphinx_dir dirname=sphinx-rootdir title="A Tutorial for the Odesolvers Package for Solving Ordinary DifferentialEquations" author="Hans Petter Langtangen and Liwei Wang" theme=$theme $name

#doconce subst '\#html_theme_options.+' "html_theme_options = {\n  'headerfont': '\"Abel\", sans-serif',\n  'bodyfont': '\"Cabin\", sans-serif',\n  'headerbg': 'black',\n  'headercolor1': 'black',}" sphinx-rootdir/conf.py

python automake-sphinx.py




