function save_figure_3x(file_prefix,varargin)

if size(varargin,1)>0
    h_fig = varargin{1};
else
    h_fig = gcf;
end
set(h_fig,'PaperPositionMode','auto');
print(h_fig,'-djpeg',file_prefix)
print(h_fig,'-depsc',file_prefix)
print(h_fig,'-dtiff',file_prefix)
print(h_fig,'-dpdf',file_prefix,'-fillpage')
saveas(h_fig,file_prefix,'fig')

end