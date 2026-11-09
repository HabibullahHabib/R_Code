from pathlib import Path
import shutil,subprocess,os
ROOT=Path(__file__).resolve().parent;FULL=ROOT/'full';OUT=FULL/'IBScovars_ReferenceStyle_Output';READY=FULL/'latex-ready';READY.mkdir(exist_ok=True)
for p in OUT.glob('*.tex'):shutil.copy2(p,READY/p.name)
shutil.copy2(FULL/'comparison_graph.pdf',READY/'IBScovars_ReferenceStyle_ForestPlot.pdf')
(READY/'Compile_Application.tex').write_text('\\documentclass{article}\n\\usepackage[a4paper,margin=1.5cm]{geometry}\n\\usepackage{graphicx}\n\\begin{document}\n\\input{IBScovars_ReferenceStyle_Manuscript_Insert.tex}\n\\end{document}\n')
with (ROOT/'latex_compile.log').open('wb') as log:r=subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error','Compile_Application.tex'],cwd=READY,stdout=log,stderr=subprocess.STDOUT)
assert r.returncode==0,(ROOT/'latex_compile.log').read_text()[-1500:]
for p in READY.iterdir():
 if p.suffix in ('.aux','.log'):p.unlink()
print('Original LaTeX exports compiled successfully.')
