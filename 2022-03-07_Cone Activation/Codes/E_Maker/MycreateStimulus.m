function E = MycreateStimulus(offset,~,~,alt,~)


% load('BAK1044R_2020_05_01_12_35_03_AOSLO.mat')
% gaussw=CFG.gaussw;
gaussw=1;
% if CFG.drift==1
%     cd([pwd,'\tempStimulus']);
% else
n_frames=1;
gap=(offset);

% % CFG.stim='subp';
% if strcmp(CFG.stim,'subp')
gap=round(gap*10);
% end

gap=gap;
padw=gap;
e=zeros(gap*5,gap*5);
e(gap+1:2*gap,gap+1:end)=1;
e(3*gap+1:4*gap,gap+1:end)=1;

%         e = imresize(e,0.01);

padh=ones(padw,size(e,2));
padv=ones(size(e,1)+2*padw,padw);

E=[padh;e;padh];
E=[padv E padv];

E=imrotate(E,alt);

% if strcmp(CFG.stim,'subp')
E = imresize(E,0.1);
f = fspecial('gaussian',10,gaussw);
E = imfilter(E,f,'replicate');
% end
canv=E;
%     end

% % CFG.flipv=0;
% if CFG.flipv == 1
%     canv=flipud(canv);
% end

% % CFG.positive=0;
% if CFG.positive == 1
%     canv=imcomplement(canv);
% end

if ~isempty(canv(canv>1)) % to avoid ripple
    canv(canv>1)=1;
end
fid = fopen('frame2.buf','w');
fwrite(fid,size(canv,2),'uint16');
fwrite(fid,size(canv,1),'uint16');
fwrite(fid, canv, 'double');
fclose(fid);
end