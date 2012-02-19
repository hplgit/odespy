#!/bin/sh

name=wrap_odesolvers

doconce format latex $name
ptex2tex -DMINTED $name
#latex -shell-escape $name
#dvipdf $name

doconce format sphinx $name
theme=pyramid
#theme=fenics_minimal
doconce sphinx_dir dirname=sphinx-rootdir title="A Tutorial for the Odesolvers Package for Solving Ordinary DifferentialEquations" author="Hans Petter Langtangen and Liwei Wang" theme=$theme $name

#doconce subst '\#html_theme_options.+' "html_theme_options = {\n  'headerfont': '\"Abel\", sans-serif',\n  'bodyfont': '\"Cabin\", sans-serif',\n  'headerbg': 'black',\n  'headercolor1': 'black',}" sphinx-rootdir/conf.py

python automake-sphinx.py




