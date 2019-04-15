fs=192000;
% fs=96000;
nfreq=30;
pathname='..\player';
figureNum=3;

global exportImage;
exportImage=0;

global hpName;

path = uigetdir(pathname);
spl=strsplit(path,'\');
folderName=char(spl(end));
hpName=folderName;
path=join([path,'\distortion\']);
mkdir(path);

answer = questdlg('l or r', ...
    '', ...
    'L','R','R');
% Handle response
switch answer
    case 'L'
        LR = 'L';
    case 'R'
        LR = 'R';
end

answer = questdlg('hd or mtd', ...
    '', ...
    'hd','mtd','hd');
% Handle response
if strcmp(answer,'hd')
    hdAnalysis(fs,nfreq,path,LR,figureNum);
elseif strcmp(answer,'mtd')
    mtdAnalysis(fs,path,LR,figureNum);
end


function mtdAnalysis(fs,path,LR,figureNum)
global exportImage;
global hpName;
sl=fs*2.5;
xl = 0;
nfreq2=30;
logfreq2=linspace(log10(20),log10(20000),nfreq2);
%     delta1k=min(abs(logfreq2-log10(1000)));
%     if min(abs(logfreq2+delta1k-log10(1000)))<delta1k
%         logfreq2=logfreq2+delta1k;
%     else
%         logfreq2=logfreq2-delta1k;
%     end
for i=logfreq2
    exfreq=10^i;
    xl=xl+sin(2*pi*exfreq/fs*(1:sl)+rand*pi)*db2mag(-30);
end
%     xl=square(2*pi*20/fs*(1:sl));
multitoneFile=strcat(path,LR,'multitone','.wav');
noisefloorFile=strcat(path,LR,'noisefloor','.wav');

[yl,fs] = audioLoader(multitoneFile,fs,xl,sl,false);
[ym,fs] = audioLoader(noisefloorFile,fs,xl*0,sl,false);
noiseLvl=mag2db(mean(abs(ym)))

figure(figureNum);
clf;
subplot(511);
plot((0:length(yl)-1)/fs,yl);
subplot(512);

if exportImage
    figure(4);
    clf;
end


[yldft, freql] = mydft(yl,fs);
yldftdb=mag2db(abs(yldft));
semilogx(freql,yldftdb);
hold on;
yln=yl;
for i=logfreq2
    exfreq=10^i;
    %         if exfreq<999 || exfreq>1001
     d = designfilt('bandstopiir','FilterOrder',2, ...
         'HalfPowerFrequency1',exfreq-i*1,'HalfPowerFrequency2',exfreq+i*1, ...
         'DesignMethod','butter','SampleRate',fs);
%d=designfilt('bandstopiir', 'FilterOrder', 2, 'StopbandFrequency1', exfreq-0.01, 'StopbandFrequency2', exfreq+0.01, 'StopbandAttenuation', 45,'SampleRate', fs);
    yln=filtfilt(d,yln);
    %         end
end

fs2=1e4;
freq2=log10(20):1/fs2:log10(2e4);

[ylndft, freql] = mydft(yln,fs);
ylndftdb=mag2db(abs(ylndft));
ylndftdb2=interp1(log10(freql(2:end)),ylndftdb(2:end),freq2);


[ymdft, freqm] = mydft(ym,fs);
ymdftdb=mag2db(abs(ymdft));


%     ymdftdb2=smooth(ymdftdb2,200)';


ptr20hz=round(20/(fs/length(yl)))+1;
ptr20000hz=round(20000/(fs/length(yl)))+1;
ylndftdb=ylndftdb(ptr20hz-1:ptr20000hz+1);
freql=freql(ptr20hz-1:ptr20000hz+1);
ymdftdb=ymdftdb(ptr20hz-1:ptr20000hz+1);
freqm=freqm(ptr20hz-1:ptr20000hz+1);

semilogx(freql,ylndftdb);
semilogx(freqm,ymdftdb);

%     ylndftdb2=smooth(ylndftdb2,10);

%     ymdftdb2=smooth(ymdftdb2,10)';

%     ylndftdb2=movmax(ylndftdb2,10);

%     ymdftdb2=movmax(ymdftdb2,10);

%     ymdftdb=smooth(ymdftdb,3)';
ymdftdb=movmax(ymdftdb,30);

ymdftdb2=interp1(log10(freqm(2:end)),ymdftdb(2:end),freq2);

mtdSpectrum=max(0,ylndftdb-ymdftdb);
%     mtdSpectrum=smooth(mtdSpectrum,2e2);
%     fs3=1e6;
%     freq3=log10(20):1/fs3:log10(2e4);
%     mtdSpectrum3=interp1(log10(freql(2:end)),mtdSpectrum(2:end),freq3);

mtdScore=mean(mtdSpectrum)
hold off;
ylim([-100 100]);
xlim([20 20000]);

if exportImage
    setFRPlot(false);
    ylabel('dB');
    if LR=='L'
        leftright='Left';
    else
        leftright='Right';
    end
    mytitle(strcat(hpName,'(',leftright,')'));
    legend('Fundamental','Noise and Distortion','Noise Floor','FontSize',15);
    savePlot(path,'mtd',LR);
    figure(figureNum);
end

subplot(513);
semilogx(freql,ylndftdb);
%     plot(freq2,ylndftdb2);
hold on;
semilogx(freqm,ymdftdb);
%     plot(freq2,ymdftdb2);
hold off;
if exportImage
    setFRPlot(false);
end
%     xlim([20 20000]);
%     ylim([-20 100]);
%     subplot(413);
%     sinad(yl,fs);
%     audioIO(yl,fs);
%     xlim([20 2e4]);
subplot(514);
semilogx(freql,(mtdSpectrum));
setFRPlot(false);
%     xlim([20 2e4]);
ylim([0 20]);
subplot(515);
%     plot(freq3,(mtdSpectrum3));
plot(freql,smooth(mtdSpectrum,2000));
ylim([0 1]);
%     xlim([log10(20) log10(2e4)]);
mtdData=[mtdScore];
csvwrite(strcat(path,'mtd',LR,'.csv'),mtdData);
end

function hdAnalysis(fs,nfreq,path,LR,figureNum)
sl=fs*1.5;
hdData=[];
for i=linspace(log10(20),log10(10000),nfreq)
    exfreq=10^i;
    xk = sin(2*pi*exfreq/fs*(1:sl))*db2mag(-0.2);
    % yj=audioIO(xk*0,fs);
    filename=strcat(path,LR,'sine',num2str(exfreq,'%.2f'),'hz.wav');
    [yk,fs]=audioLoader(filename,fs,xk,sl,true);
    figure(figureNum);
    clf;
    subplot(411);
    plot((0:length(yk)-1)/fs,yk);
    subplot(412);
    % yj=yj(0.5*fs:sl);
    %         yk=yk(0.5*fs:sl);
    
    %         audiowrite(strcat(path,LR,'sine',num2str(exfreq,'%.2f'),'hz.wav'),yk,fs);
    
    % [yjdft, freqj] = mydft(yj,fs);
    [ykdft, freqk] = mydft(yk,fs);
    semilogx(freqk,mag2db(abs(ykdft)));
    % hold on;
    % semilogx(freqj,mag2db(abs(yjdft)));
    % hold off;
    xlim([20 20000]);
    subplot(413);
    thd(yk,fs,min(40,floor(20000/exfreq)));
    xlim([20/1000 20000/1000]);
    
    [r,harmpow,harmfreq]= thd(yk,fs,min(40,floor(20000/exfreq)));
    
    thdArray=[exfreq r -snr(yk,fs) harmpow'];
    if isempty(hdData)
        hdData=thdArray;
    else
        if length(thdArray)<length(hdData(1,:))
            thdArray(length(hdData(1,:)))=0;
        end
        hdData=[hdData;thdArray];
    end
    subplot(414);
    snr(yk,fs);
    xlim([20/1000 20000/1000]);
end
csvwrite(strcat(path,'hd',LR,'.csv'),hdData);
end

function [yl,fsnew]=audioLoader(filename,fs,xl,sl,isSine)
if isfile(filename)
    [yl,fsnew] = audioread(filename);
    assert(fs==fsnew);
else
    yl= audioIO(xl,fs);
    yl=yl(0.5*fs:sl);
    audiowrite(filename,yl,fs,'BitsPerSample',24);
    
    fsnew=fs;
end
end
