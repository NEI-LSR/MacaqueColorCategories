# MacaqueColorCategories

This repo is the home for the Macaque Color Categories project.

## Manuscript

Overleaf link: https://www.overleaf.com/9918797415dwvhydynjtmv (it's private so you'll need to be invited)

To convert the LaTeX files to word, use pandoc (from the writing/LaTeX folder):
```
pandoc -s main.tex --citeproc --bibliography=bib.bib -o ../pandocOutput.docx
```
and then rename and put the old version in `old`.

## Code requirements

There are various submodules on this project.
To fill them, either clone the repo with `git clone --recurse-submodules -j8 git://github.com/foo/bar.git` or clone the repo normally and then run:
```
git submodule init
git submodule update
```

This project also requires PsychToolbox (though not super heavily, mainly just LuvToXYZ to plotting, I think...)
