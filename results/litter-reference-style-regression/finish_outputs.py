from pathlib import Path
import os,subprocess,shutil,csv
ROOT=Path(__file__).resolve().parent;FULL=ROOT/'full'
env_repair=os.environ.copy();env_repair.update(PATH='/workspace/.r-runtime/bin:'+env_repair['PATH'],R_LIBS_SITE='/workspace/.r-runtime/usr/lib/R/site-library',R_LIBS_USER='/workspace/.r-runtime/library',LD_LIBRARY_PATH='/workspace/.r-runtime/usr/lib/x86_64-linux-gnu')
subprocess.run(['/workspace/.r-runtime/bin/Rscript','--vanilla',str(ROOT/'repair_validation_summary.R'),str(ROOT)],env=env_repair,check=True)
subprocess.run(['python',str(ROOT/'make_report.py')],check=True)
env=os.environ.copy();env.update(PATH='/workspace/.r-runtime/bin:'+env['PATH'],R_LIBS_SITE='/workspace/.r-runtime/usr/lib/R/site-library',R_LIBS_USER='/workspace/.r-runtime/library',LD_LIBRARY_PATH='/workspace/.r-runtime/usr/lib/x86_64-linux-gnu',XDG_CACHE_HOME='/workspace/.r-runtime/cache')
p=subprocess.run(['/workspace/.r-runtime/bin/Rscript','--vanilla',str(ROOT/'render_highlighted_plot.R'),str(ROOT)],env=env,capture_output=True,text=True);assert p.returncode==0,p.stderr
subprocess.run(['pdfunite',str(FULL/'main_tables.pdf'),str(FULL/'comparison_graph.pdf'),str(FULL/'simulation_report.pdf')],check=True)
rows=list(csv.DictReader((FULL/'Litter_Application_Final/litter_PRIMARY_REG_APPLICATION_TABLE.csv').open()));ready=FULL/'latex-ready';ready.mkdir(exist_ok=True)
tex=[r'\begin{table}[htbp]',r'\centering',r'\caption{Litter application: regression imputation, 40 imputations.}',r'\begin{tabular}{lccc}',r'\hline',r'Dose comparison & Point estimate & Simultaneous 95\% CI & Length \\',r'\hline']
for method in ['Proposed method with shrinkage','Proposed method without shrinkage','Dunnett-type procedure with naive case deletion']:
 tex.append(r'\multicolumn{4}{l}{'+method+r'} \\')
 for x in rows:
  if x['method']==method:tex.append(f"{x['contrast']} & {float(x['estimate']):.3f} & ${x['CI']}$ & {float(x['length']):.3f} "+r'\\')
tex += [r'\hline',r'\end{tabular}',r'\end{table}'];(ready/'Litter_Regression_Table.tex').write_text('\n'.join(tex)+'\n')
shutil.copy2(FULL/'comparison_graph.pdf',ready/'comparison_graph.pdf')
figure='\n'.join([r'\begin{figure}[p]',r'\centering',r'\includegraphics[width=0.98\textwidth]{comparison_graph.pdf}',r'\caption{Regression-imputation application. Highlighted rows show the proposed with-shrinkage and without-shrinkage results combining the same 40 imputations using Rubin rules.}',r'\end{figure}']);(ready/'Litter_Regression_Figure.tex').write_text(figure+'\n')
(ready/'Compile_Application.tex').write_text('\n'.join([r'\documentclass{article}',r'\usepackage[a4paper,margin=1.5cm]{geometry}',r'\usepackage{graphicx}',r'\begin{document}',r'\input{Litter_Regression_Table.tex}',r'\input{Litter_Regression_Figure.tex}',r'\end{document}'])+'\n')
with (ROOT/'latex_compile.log').open('wb') as log:p=subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error','Compile_Application.tex'],cwd=ready,stdout=log,stderr=subprocess.STDOUT)
assert p.returncode==0,(ROOT/'latex_compile.log').read_text()[-1500:]
for p in ready.iterdir():
 if p.suffix in ('.aux','.log'):p.unlink()
print('PDF, highlighted graph and compiled LaTeX completed.')
