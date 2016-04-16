#!/bin/sh
set -x

sh clean.sh
doconce spellcheck -d .dict4spell.txt *.do.txt
if [ $? -ne 0 ]; then
   echo "Spelling errors - abort."
   exit 1
fi

name=main_odespy

doconce format pdflatex $name "--latex_code_style=default:lst[style=yellow2_fb]@sys:lst-gray"
pdflatex $name
pdflatex $name

doconce format html $name --html_style=bootswatch_readable --html_code_style=inherit --html_output=odespy
doconce format html $name --html_style=solarized3 --html_output=odespy-solarized

doconce format sphinx $name
theme=cbc
#theme=fenics_minimal
doconce sphinx_dir dirname=sphinx-rootdir title="A Tutorial for the Odespy Interface to ODE Solvers" copyright="Hans Petter Langtangen and Liwei Wang" theme=$theme $name

#doconce subst '\#html_theme_options.+' "html_theme_options = {\n  'headerfont': '\"Abel\", sans-serif',\n  'bodyfont': '\"Cabin\", sans-serif',\n  'headerbg': 'black',\n  'headercolor1': 'black',}" sphinx-rootdir/conf.py

python automake_sphinx.py

# Publish
dest=../../pub/tutorial
cp -r sphinx-rootdir/_build/html $dest
cp $name.pdf $dest/odespy.pdf
cp -r fig-odespy $dest
cp odespy*.html $dest

# Add to git if new files have been created
git add $dest
