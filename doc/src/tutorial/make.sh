#!/bin/sh

name=odesolvers_wrap

doconce format latex $name
ptex2tex -DMINTED $name
#latex -shell-escape $name
#dvipdf $name

doconce format sphinx $name
theme=agogo
theme=fenics_minimal
doconce sphinx_dir dirname=sphinx-rootdir theme=$theme $name
cp -r ex sphinx-rootdir  # figures
#doconce subst '\#html_theme_options.+' "html_theme_options = {\n  'headerfont': '\"Abel\", sans-serif',\n  'bodyfont': '\"Cabin\", sans-serif',\n  'headerbg': 'black',\n  'headercolor1': 'black',}" sphinx-rootdir/conf.py

sh automake-sphinx.sh



