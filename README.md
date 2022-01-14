# MacaqueColorCategories

This repo is the home for a manuscript on the topic of color categories in macaque.

Overleaf link: https://www.overleaf.com/9918797415dwvhydynjtmv (it's private so you'll need to be invited)

Don't delete the .git hidden folder.

To convert the LaTeX files to word, use pandoc (from the LaTeX folder):
```
pandoc main.tex --bibliography=bib.bib -o ../pandocOutput.docx
```
and then rename and put the old version in `old`.