function res = save2pdf(name,varargin)
pathToPdflatex = '/Library/TeX/texbin/pdflatex' ; 

files = cell(size(varargin)) ;
for ii = 1:numel(varargin)
    files{ii} = sprintf('%s_fig%g.pdf',name,ii) ;
    print(varargin{ii},'-dpdf','-painters',files{ii}) ;
end

fh = fopen(sprintf('%s.tex',name),'w+') ;
fprintf(fh,'\\documentclass{article}\n') ;
fprintf(fh,'\\usepackage{graphicx}\n') ;
fprintf(fh,'\\begin{document}\n') ;
for ii = 1:numel(files)
    fprintf(fh,'\\includegraphics[width=\\textwidth]{%s}\n\\newpage\n',files{ii}) ;
end
fprintf(fh,'\\end{document}\n') ;
fclose(fh) ;

[~,res] = system(sprintf('%s %s.tex',pathToPdflatex,name)) ;
disp(res)
end