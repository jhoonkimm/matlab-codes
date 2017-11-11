function [mov, loadedAVI] = playAVI(aviFilename)
%%

loadedAVI = VideoReader(aviFilename);

mov(loadedAVI.NumberOfFrames) = struct('cdata',[],'colormap',[]);
for ii = 1:loadedAVI.NumberOfFrames
    mov(ii) = im2frame(read(loadedAVI,ii));
end

set(gcf,'position', [150 150 loadedAVI.Width loadedAVI.Height])
set(gca,'units','pixels');
set(gca,'position',[0 0 loadedAVI.Width loadedAVI.Height])

image(mov(1).cdata,'Parent',gca);
axis off;

movie(mov, 1, loadedAVI.FrameRate);

end

