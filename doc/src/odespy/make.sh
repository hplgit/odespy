#!/bin/sh
doconce spellcheck -d .dict4spell.txt *.do.txt
if [ $? -ne 0 ]; then
   echo "Spelling errors - abort."
   exit 1
fi

name=main_odespy

doconce format pdflatex $name
ptex2tex -DMINTED -DLATEX_HEADING=traditional $name
pdflatex -shell-escape $name
pdflatex -shell-escape $name

doconce format sphinx $name
theme=pyramid
#theme=fenics_minimal
doconce sphinx_dir dirname=sphinx-rootdir title="A Tutorial for the Odespy Interface to ODE Solvers" author="Hans Petter Langtangen and Liwei Wang" theme=$theme $name

#doconce subst '\#html_theme_options.+' "html_theme_options = {\n  'headerfont': '\"Abel\", sans-serif',\n  'bodyfont': '\"Cabin\", sans-serif',\n  'headerbg': 'black',\n  'headercolor1': 'black',}" sphinx-rootdir/conf.py

python automake_sphinx.py

doconce format html $name

# Copy to tutorial
cp -r sphinx-rootdir/_build/html ../../tutorial/
cp $name.pdf ../../tutorial/odespy.pdf
cp -r fig-odespy ../../tutorial/
cp $name.html ../../tutorial/odespy.html
doconce format html $name --html_style=solarized
cp $name.html ../../tutorial/odespy-solarized.html

