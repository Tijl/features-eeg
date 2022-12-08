%%


%%
f = figure(1);clf
f.Position(3:4) = [1000 400];
% stimuli

im = [];
fns = dir('../stimuli_small/*.png');
assert(numel(fns)==256);
for f=1:256
    fn = fullfile(fns(f).folder,fns(f).name);
    [i,map] = imread(fn);
    if ndims(i)<3 %(some ims were saved as indexed, so make them all truecolour)
        i = ind2rgb(i,map);
    else
        i = im2double(i);
    end
    im(:,:,:,f) = i;
end

a=subplot(1,2,1)
for i=1:256
    aw = .5;
    x = floor((i-1)/16);
    y = mod((i-1),16);
    c = im(:,:,:,i);
    image('XData',x+[-aw aw],'YData',y+[-aw aw],'CData',c)
end
xlim([-.5 16.5])
ylim([-.5 16.5])
axis off
%montage(im)
for i=1:4:13
    text(-1,i,'<- SF ->','HorizontalAlignment','center','Rotation',90,'FontSize',14)
    text(i,-1,'<- Colour ->','HorizontalAlignment','center','FontSize',14)
end
text(7,-2,'<- Contrast ->','HorizontalAlignment','center','FontSize',15)
text(-2,7,'<- Orientation ->','HorizontalAlignment','center','Rotation',90,'FontSize',15)

a=subplot(2,4,3)
rng(1);
idx = randsample(1:256,10);
for i=10:-1:1
    aw = .5;
    x = i*.6;
    y=i*.4;
    c = im(:,:,:,idx(i));
    if ~mod(i,2)
        c(:) = c(1);
    end
    image('XData',x+[-aw aw],'YData',y+[-aw aw],'CData',c)
end
axis equal
axis off


a=subplot(2,4,7)
idx = randsample(1:256,3);
for i=1:3
    aw = .5;
    x = .6*i;
    y=(i==2)*1.2;
    image('XData',x+[-aw aw],'YData',y+[-aw aw],'CData',im(:,:,:,idx(i)))
end
axis equal
axis off


a=annotation('textbox',[.1,.95,.2,.01],'string','A: Stimuli','HorizontalAlignment','left','LineStyle','none','FontWeight','bold','FontSize',15)
a=annotation('textbox',[.5,.95,.2,.01],'string','B: EEG Design','HorizontalAlignment','left','LineStyle','none','FontWeight','bold','FontSize',15)
a=annotation('textbox',[.5,.55,.2,.01],'string','C: online triplet task','HorizontalAlignment','left','LineStyle','none','FontWeight','bold','FontSize',15)

%% save
fn = '../figures/design';
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im=imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');


