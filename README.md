# MacaqueColorCategories

This repo is the home for the Macaque Color Categories project.

Overleaf link: https://www.overleaf.com/9918797415dwvhydynjtmv (it's private so you'll need to be invited)

Don't delete the `.git` hidden folder.

To convert the LaTeX files to word, use pandoc (from the writing/LaTeX folder):
```
pandoc -s main.tex --citeproc --bibliography=bib.bib -o ../pandocOutput.docx
```
and then rename and put the old version in `old`.
