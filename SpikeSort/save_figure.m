function save_figure(filename,directory,suffix,varargin)

% parse inputs
p = inputParser;
default_jpg = 0;
default_fig = 0;
default_tif = 0;
default_pdf = 0;
default_eps = 0;
default_svg = 0;
default_h_fig = gcf;

addOptional(p,'jpg',default_jpg,@isnumeric);
addOptional(p,'fig',default_fig,@isnumeric);
addOptional(p,'tif',default_tif,@isnumeric);
addOptional(p,'pdf',default_pdf,@isnumeric);
addOptional(p,'eps',default_eps,@isnumeric);
addOptional(p,'svg',default_svg,@isnumeric);
addOptional(p,'h_fig',default_h_fig);

parse(p,varargin{:});

% parse path stuff
[pathstr,record_name,ext] = fileparts(filename);
% check if intended Figures subdirectory exists, if not create it
save_path = fullfile(pathstr,'Figures',directory);
if ~exist(save_path,'file')
    mkdir(save_path)
end
save_filename = fullfile(save_path,strcat(record_name,'_',suffix));

% save stuff!
%set(p.Results.h_fig,'PaperPositionMode','auto');
if p.Results.jpg
    print(p.Results.h_fig,'-djpeg',save_filename)
end
if p.Results.eps
    print(p.Results.h_fig,'-depsc',save_filename)
end
if p.Results.tif
    print(p.Results.h_fig,'-dtiff',save_filename)
end
if p.Results.pdf
    print(p.Results.h_fig,'-dpdf',save_filename)
end
if p.Results.svg
    saveas(p.Results.h_fig, [save_filename '.svg'])
end
if p.Results.fig
    saveas(p.Results.h_fig,save_filename,'fig')
end
end