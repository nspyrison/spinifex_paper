
tools::texi2pdf(file="paper/RJwrapper.tex",clean = T)
## Error in texi2dvi(file = file, pdf = TRUE, clean = clean, quiet = quiet,  :
## unable to run 'pdflatex' on 'paper/spinifex-paper.tex'

Sys.which("pdflatex") # tinytex...pdflatex.exe
Sys.getenv("PATH") #. large list

### I think this path is not the same use case.
library(msa)
#msa::msaPrettyPrint("paper/spinifex-paper.tex", output="tex", showConsensus = "none", askForOverwrite=FALSE, verbose=FALSE,
#               file = "tex_file")
?msa::msaPrettyPrint

## read sequences
filepath <- system.file("examples", "exampleAA.fasta", package="msa")
mySeqs <- readAAStringSet(filepath)

## call unified interface msa() for default method (ClustalW) and
## default parameters
myAlignment <- msa(mySeqs)

## show resulting LaTeX code with default settings
msaPrettyPrint(myAlignment, output="asis", askForOverwrite=FALSE)

## create PDF file according to some custom settings
tmpFile <- tempfile(pattern="msa", tmpdir=".", fileext=".pdf")
tmpFile
msaPrettyPrint(myAlignment, file=tmpFile, output="pdf",
               showNames="left", showNumbering="none", showLogo="top",
               showConsensus="bottom", logoColors="rasmol",
               verbose=FALSE, askForOverwrite=FALSE)

