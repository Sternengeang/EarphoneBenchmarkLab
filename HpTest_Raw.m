for i=1:4
    figure(i);
    clf;
end

global exportImage;
exportImage=0;
global showFigure;
showFigure=0;
global isPhaseFlipped;
isPhaseFlipped=0;
global isPhase360;
isPhase360=[0,0];

runFolder=0;

resultFile='result_2.0raw.csv';
result=readtable(resultFile);
if runFolder
    D = dir('..\5_2\_chosen');
    for k = 3:length(D)
        if D(k).isdir
            [wavL,wavR,hdL,hdR,mtdL,mtdR,folderName,path]=readWav(strcat(D(k).folder,'\',D(k).name));
            result=folderAnalysis(wavL,wavR,hdL,hdR,mtdL,mtdR,folderName,path,result);
        end
    end
else
    [wavL,wavR,hdL,hdR,mtdL,mtdR,folderName,path]=readWav('');
    result=folderAnalysis(wavL,wavR,hdL,hdR,mtdL,mtdR,folderName,path,result);
end

writetable(result,resultFile);

function result=folderAnalysis(wavL,wavR,hdL,hdR,mtdL,mtdR,folderName,path,result)
global exportImage
global freq2;
global hpName;
global ptr10000hz;
global isPhaseFlipped;
global curve2;
global isPhase360;
if ~any(ismember(result.name,folderName))
    result=[result;{folderName,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}];
end
[fr2L,ph2L,result]=singleChannelAnalysis(wavL,hdL,mtdL,1,folderName,path,result);
[fr2R,ph2R,result]=singleChannelAnalysis(wavR,hdR,mtdR,2,folderName,path,result);
lrFrDiff=0;
for i=1:length(fr2R)
    if i<=ptr10000hz
        k=2;
    else
        k=1;
    end
    lrFrDiff=lrFrDiff+k*abs(fr2L(i)-fr2R(i));
end
lrFrDiff=(lrFrDiff/(ptr10000hz*2+length(fr2R)-ptr10000hz));
lrFrDiff=-18*lrFrDiff+100;

result(ismember(result.name,{folderName}),:).LRFRCorr=lrFrDiff;

lrPhDiff=getPhaseDiff(ph2L,ph2R,false,2);
lrPhDiff=-100/pi*lrPhDiff+100;
if(lrPhDiff>620/7)
    lrPhDiff=7*lrPhDiff-600;
else
    lrPhDiff=20/(620/7)*lrPhDiff;
end
result(ismember(result.name,{folderName}),:).LRPhaseCorr=lrPhDiff;

if exportImage
    figure(4);
    clf;
    plot(freq2,fr2L,'LineWidth',1.5);
    hold on;
    plot(freq2,fr2R,'LineWidth',1.5);
    plot(freq2,curve2,'k--','LineWidth',1.5);
    hold off;
    mytitle(strcat(hpName,''));
    setFRPlot(true);
    ylim([-50 0]);
    ylabel('dB');
    legend('Left','Right','Target','FontSize',15);
    savePlot(path,'frdiff','');
    figure(3);
    clf;
    if isPhase360(1)
        ph2L=ph2L+2*pi;
    end
    if isPhase360(2)
        ph2R=ph2R+2*pi;
    end
    if ~isPhaseFlipped
        plot(freq2,wrapToPi(ph2L)/pi*180,'LineWidth',1.5);
    else
        plot(freq2,(ph2L)/pi*180,'LineWidth',1.5);
    end
    hold on;
    if ~isPhaseFlipped
        plot(freq2,wrapToPi(ph2R)/pi*180,'LineWidth',1.5);
    else
        plot(freq2,(ph2R)/pi*180,'LineWidth',1.5);
    end
    hold off;
    mytitle(strcat(hpName,''));
    setFRPlot(true);
    ylim([-360 360]);
    ylabel('Degree');
    legend('Left Phase','Right Phase','FontSize',15);
    savePlot(path,'phdiff','');
end

frJitter=(result(ismember(result.name,{folderName}),:).SmoothL+result(ismember(result.name,{folderName}),:).SmoothR)/2;
frJitter=-20*frJitter+100;
if frJitter>50
    frJitter=8/5*frJitter-60;
else
    frJitter=20/50*frJitter;
end

phaseDiff=(result(ismember(result.name,{folderName}),:).MinPhaseL+result(ismember(result.name,{folderName}),:).MinPhaseR)/2;
phaseDiff=-100/pi*phaseDiff+100;
if phaseDiff>90
    phaseDiff=8*phaseDiff-700;
else
    phaseDiff=20/90*phaseDiff;
end

thdScore=(result(ismember(result.name,{folderName}),:).THDL+result(ismember(result.name,{folderName}),:).THDR)/2;
thdScore=100-thdScore*100;
if thdScore>99.52
    thdScore=500/3*thdScore-49700/3;
else
    thdScore=20/99.52*thdScore;
end

mtdScore=(result(ismember(result.name,{folderName}),:).MTDL+result(ismember(result.name,{folderName}),:).MTDR)/2;
mtdScore=-10*mtdScore+100;
if mtdScore>80
    mtdScore=4.5*mtdScore-350;
else
    mtdScore=10/80*mtdScore;
end

balanceScore=(result(ismember(result.name,{folderName}),:).FRBalanceL+result(ismember(result.name,{folderName}),:).FRBalanceR)/2;
balanceScore=-100/40*balanceScore+100;

referenceDiff=(result(ismember(result.name,{folderName}),:).ReferenceL+result(ismember(result.name,{folderName}),:).ReferenceR)/2;
referenceDiff=-10*referenceDiff+100;

result(ismember(result.name,{folderName}),:).SmoothAvg=frJitter;
result(ismember(result.name,{folderName}),:).MinPhaseAvg=phaseDiff;
result(ismember(result.name,{folderName}),:).THDAvg=thdScore;
result(ismember(result.name,{folderName}),:).MTDAvg=mtdScore;
result(ismember(result.name,{folderName}),:).FRBalanceAvg=balanceScore;
result(ismember(result.name,{folderName}),:).ReferenceAvg=referenceDiff;


frFinalScore=(result(ismember(result.name,{folderName}),:).SmoothAvg*1 ...
    +result(ismember(result.name,{folderName}),:).ReferenceAvg*4 ...
    +result(ismember(result.name,{folderName}),:).LRFRCorr*1 ...
    )/(1+4+1);
phFinalScore=(result(ismember(result.name,{folderName}),:).MinPhaseAvg*1 ...
    + result(ismember(result.name,{folderName}),:).LRPhaseCorr*3 ...
    )/4;
hdFinalScore=result(ismember(result.name,{folderName}),:).THDAvg;
mtdFinalScore=result(ismember(result.name,{folderName}),:).MTDAvg;
stereoFinalScore=(result(ismember(result.name,{folderName}),:).LRFRCorr+result(ismember(result.name,{folderName}),:).LRPhaseCorr)/2;
result(ismember(result.name,{folderName}),:).FRAvg=frFinalScore;
result(ismember(result.name,{folderName}),:).StereoAvg=stereoFinalScore;
result(ismember(result.name,{folderName}),:).PHAvg=phFinalScore;
total=(frFinalScore*25+phFinalScore*2+hdFinalScore*3+mtdFinalScore*1)/(25+2+3+1);
if exportImage
    scores=[frFinalScore phFinalScore hdFinalScore mtdFinalScore];
    weight=[65 10 15 5  ];
    myPolarPlot(scores,weight,total,folderName,path);
end
result(ismember(result.name,{folderName}),:).Total=total;
end

function [fr2,phase2,result]=singleChannelAnalysis(wavFile, hdFile,mtdFile, figureNum,folderName,path,result)
global exportImage;
global hpName;
global isPhaseFlipped;
global showFigure;
global curve2;
global freq2;
global ptr10000hz;
global isPhase360;
fs2=1e2;
fr2=[];
phase2=[];
if figureNum==1
    LR='L';
    leftright='Left';
elseif figureNum==2
    LR='R';
    leftright='Right';
end
figure(figureNum);
if isfile(wavFile)
    [y,fs] = audioread(wavFile);
    curve=csvread('new target.csv');
    
    ydft = fft(y);
    subplot(527);
    yi=ifft(ydft);
    if showFigure
        plot((1:length(yi))/fs,yi);
        grid;
    end
    ydft = ydft(1:length(y)/2+1);
    freq = 0:fs/length(y):fs/2;
    subplot(521);

    freq2=log10(20):1/fs2:log10(2e4);
    fr2=interp1(log10(freq(2:end)),mag2db(abs(ydft(2:end))),freq2);
    if showFigure
        plot(freq2,fr2);
        hold on;
    end
    curveAlignFreq=1000;
    distTo1k=10000;
    ptrTo1k=0;
    for i=1:length(curve(:,1))
        tmpDist=abs(curve(i,1)-curveAlignFreq);
        if tmpDist<distTo1k
            distTo1k=tmpDist;
            ptrTo1k=i;
        end
    end
    curve1000=curve(ptrTo1k,2);
    for i=1:length(freq)-1
        if freq(i)==curveAlignFreq
            ydft1000=mag2db(abs(ydft(i)));
        elseif freq(i)<curveAlignFreq && freq(i+1)>curveAlignFreq
            ydft1000=mag2db((abs(ydft(i))+abs(ydft(i+1)))/2);
        end
    end
    curve2=interp1(log10(curve(:,1)),curve(:,2),freq2)-curve1000+ydft1000;
    if showFigure
        plot(freq2,curve2);
        hold off;
        axis([log10(20) log10(20000) -60 0]);
        grid;
        title('Frequency Response');
        xlabel('Hz');
        ylabel('dB');
    end
    
    
    if exportImage
        figure(3);
        plot(freq2,fr2,'LineWidth',1.5);
        hold on;
        plot(freq2,curve2,'k--','LineWidth',1.5);
        hold off;
        mytitle(strcat(hpName,''));
        setFRPlot(true);
        ylim([-50 0]);
        ylabel('dB');
        legend(leftright,'Target','FontSize',15);
        savePlot(path,'fr',LR);
        figure(4);
        plot(freq2,fr2-curve2-curve1000+ydft1000,'LineWidth',1.5);
        mytitle(strcat(hpName,'(In-Ear Target)'));
        legend(leftright,'FontSize',15);
        setFRPlot(true);
        ylim([-50 0]);
        ylabel('dB');
        savePlot(path,'fr-curve',LR);
        figure(figureNum);
    end
    
    ptr50hz=round((log10(50)-log10(20))*fs2);
    ptr100hz=round((2-log10(20))*fs2);
    ptr1000hz=round((3-log10(20))*fs2);
    
    ptr10000hz=round((4-log10(20))*fs2);
    ptr15000hz=round((log10(15000)-log10(20))*fs2);
    ptr16000hz=round((log10(16000)-log10(20))*fs2);
    ptr200hz=fs2ptr(fs2,200);
    ptr8000hz=fs2ptr(fs2,8000);
    ptr20000hz=fs2ptr(fs2,20000);
    
    db50hz=fr2(ptr50hz);
    db100hz=fr2(ptr100hz);
    
    
    
    frJitter=0;
    fr2comp=fr2(1:ptr20000hz)-curve2(1:ptr20000hz);
    fr2smooth=smooth(fr2comp(1:ptr16000hz),35);
    
    for i=1:ptr1000hz
        frJitter=frJitter+abs(fr2comp(i)-fr2smooth(i))*1;
    end

    for i=ptr1000hz+1:ptr8000hz
        frJitter=frJitter+abs(fr2comp(i)-fr2smooth(i))*4;
    end
    for i=ptr8000hz+1:ptr16000hz
        frJitter=frJitter+abs(fr2comp(i)-fr2smooth(i))*1;
    end
    frJitter=frJitter/(ptr1000hz*1 ...
        +(ptr8000hz-ptr1000hz)*4 ...
        +(ptr16000hz-ptr8000hz)*1);
    
    if figureNum==1
        result(ismember(result.name,{folderName}),:).SmoothL=frJitter;
    else
        result(ismember(result.name,{folderName}),:).SmoothR=frJitter;
    end    
    
    fr2Groups={
        1:ptr50hz;
        ptr50hz+1:ptr200hz;
        ptr200hz+1:ptr1000hz;
        ptr1000hz+1:ptr8000hz;
        ptr8000hz+1:ptr16000hz;
        ptr16000hz+1:ptr20000hz;
        };
    fr2GroupWeights=[4 8 10 15 8 0.1];
    for i=1:length(fr2Groups)
        fr2compAvgs(i)=mean(fr2comp(fr2Groups{i}));
    end
    
    minBalanceScore=65535;
    minMean=65535;
    minScores=[];
    for imean=min(fr2comp):0.1:max(fr2comp)
        balanceDeltas=(fr2compAvgs-imean);
        scores=[];
        balanceScore=0;
        balanceWeight=0;
        for i=1:length(balanceDeltas)
            scores(i)=abs(balanceDeltas(i))*length(fr2Groups{i})*fr2GroupWeights(i);
            balanceScore=balanceScore+abs(balanceDeltas(i))^2*length(fr2Groups{i})*fr2GroupWeights(i);
            balanceWeight=balanceWeight+length(fr2Groups{i})*fr2GroupWeights(i);
        end
        balanceScore=balanceScore/balanceWeight;
        if balanceScore<minBalanceScore
            minBalanceScore=balanceScore;
            minMean=imean;
            minBalanceDeltas=balanceDeltas;
            minScores=scores;
        end
    end
    
    
    if figureNum==1
        result(ismember(result.name,{folderName}),:).FRBalanceL=minBalanceScore;
    else
        result(ismember(result.name,{folderName}),:).FRBalanceR=minBalanceScore;
    end
    
    minReferenceDiff=0;
    referenceWeight=0;
    minMean=0;
    for i=1:length(fr2Groups)
        minScores(i)=sum(abs(fr2comp(fr2Groups{i})))*fr2GroupWeights(i);
        minReferenceDiff=minReferenceDiff+sum(abs(fr2comp(fr2Groups{i})))*fr2GroupWeights(i);
        referenceWeight=referenceWeight+length(fr2Groups{i})*fr2GroupWeights(i);
    end
    minReferenceDiff=minReferenceDiff/referenceWeight;
    
    
    if figureNum==1
        result(ismember(result.name,{folderName}),:).ReferenceL=minReferenceDiff;
    else
        result(ismember(result.name,{folderName}),:).ReferenceR=minReferenceDiff;
    end
    
    
    subplot(522);
    phase=unwrap(angle(ydft));
    [deltatime,isFlipped]=getPeakDelta(y,fs);
    phase=phaseOffset(freq,phase,deltatime);
    
    if showFigure
        semilogx(freq,(phase)/pi*180);
        xlim([20 20000]);
        ylim([-180 180]);
        title('Phase Response');
        xlabel('Hz');
        ylabel('Degree');
        hold on;
    end
    
    [c, ym]=rceps(y);
    ymdft = fft(ym);
    ymdft = ymdft(1:length(ym)/2+1);
    freqm = 0:fs/length(ym):fs/2;
    
    phasem=unwrap(angle(ymdft));
    deltatimem=getPeakDelta(ym,fs);
    phasem=phaseOffset(freqm,phasem,deltatimem);
    
    if showFigure
        semilogx(freqm,(phasem)/pi*180);
        xlim([20 20000]);
        grid;
        hold off;
    end
    
    if exportImage
        figure(3);
        if ~isPhaseFlipped
            phdisp=(phase)/pi*180;
        else
            phdisp=unwrap(phase)/pi*180;
        end
        if isPhase360(figureNum)
            phdisp=phdisp+360;
        end
        if LR=='L'
            myColor=[         0    0.4470    0.7410];
        else
            myColor=[    0.8500    0.3250    0.0980];
        end
        semilogx(freq,phdisp,'LineWidth',1.5,'Color',myColor);
        hold on;
        semilogx(freqm,(phasem)/pi*180,'k','LineWidth',1.5);
        hold off;
        setFRPlot(false);
        ylim([-360 360]);
        ylabel('Degree');
        mytitle(hpName);
        legend(strcat(leftright,' Phase'),'Minimum Phase','FontSize',15);
        savePlot(path,'phase',LR);
        figure(figureNum);
    end
    
    phase=unwrap(phase);
    phase2=interp1(log10(freq(2:end)),phase(2:end),freq2);
    phasem2=interp1(log10(freq(2:end)),phasem(2:end),freq2);
    
    
    phaseDiff=getPhaseDiff(phase2,phasem2,isFlipped,1);
    
    if figureNum==1
        result(ismember(result.name,{folderName}),:).MinPhaseL=phaseDiff;
    else
        result(ismember(result.name,{folderName}),:).MinPhaseR=phaseDiff;
    end
    
    subplot(523);
    
    if showFigure
        plot((0:length(y)-1)/fs+1+deltatime,y);
        hold on;
        plot((0:length(ym)-1)/fs+1+deltatimem,ym);
        hold off;
        xlim([1-0.001,1.005]);
        grid;
        title('Impulse Response');
        xlabel('sec');
        ylabel('% FS');
    end
    
    subplot(524);
    
    if showFigure
        gd=groupDelay(freq,phase);
        semilogx(freq,gd);
        hold on;
        gdm=groupDelay(freqm,phasem);
        semilogx(freqm,gdm);
        hold off;
        xlim([20 20000]);
        ylim([-0.02 0.02]);
        grid;
        title('Group Delay');
        xlabel('Hz');
        ylabel('sec');
    end
    
    subplot(525);
    pd=phaseDelay(freq,phase);
    
    if showFigure
        semilogx(freq,pd);
        hold on;
        pdm=phaseDelay(freqm,phasem);
        semilogx(freqm,pdm);
        hold off;
        xlim([20 20000]);
        grid;
        title('Phase Delay');
        xlabel('Hz');
        ylabel('sec');
    end
    
    if showFigure
        subplot(526);
        mycsd(y,fs,0.001/2,0.5,-deltatime,-deltatime+0.005);
        subplot(528);
        mycsd(ym,fs,0.001/2,0.5,max([1/fs -deltatimem]),max([0.005 -deltatimem+0.005]));
    end
    
end

subplot(529);
if isfile(hdFile)
    hd=csvread(hdFile);
    
    if showFigure
        semilogx(hd(:,1),hd(:,2));
        hold on;
        semilogx(hd(:,1),hd(:,3));
        hold off;
        xlim([20 20000]);
        grid;
        legend('thd','snr');
    end
    
    
    if exportImage
        figure(3);
        semilogx(hd(:,1),-hd(:,3),'LineWidth',1.5);
        setFRPlot(false);
        ylabel('dB');
        ylim([0 120]);
        legend('SNR','FontSize',15);
        mytitle(hpName);
        savePlot(path,'snr',LR);
        figure(figureNum);
    end
    
    orders=2:40;
    thd3=[];
    thd3disp=[];
    dbWeight=2;
    for i=orders
        hdsum=0;
        hdcnt=0;
        
        hddispsum=0;
        
        octavesum=0;
        octavecnt=0;
        
        j=i+3;
        for k=1:length(hd(:,1))
            if hd(k,j)~=0
                deltaDB=dbWeight*min((i-2),20);
                hdsum=hdsum+db2mag(hd(k,j)-hd(k,4)+deltaDB);
                hddispsum=hddispsum+db2mag(hd(k,j)-hd(k,4));
                hdcnt=hdcnt+1;
            end
        end
        if hdsum~=0
            hdsum=hdsum/hdcnt;
            if 2^nextpow2(i)~=i
                thd3=[thd3 hdsum];
            end
        end
        
        if hddispsum~=0
            hddispsum=hddispsum/hdcnt;
            thd3disp=[thd3disp hddispsum];
        end
        
    end
    
    if exportImage
        figure(4);
        b=bar(orders,mag2db(thd3disp),'BaseValue',-120,'FaceColor','flat');
        hold on;
        bar([2],[-120],'BaseValue',-120,'FaceColor','flat');
        hold off;
        for i=orders
            if 2^nextpow2(i)~=i
                b.CData(i-1,:)=[0.8500    0.3250    0.0980];
            end
        end
        xticks(orders);
        yticks(-120:20:0);
        ylim([-120,0]);
        mytitle(strcat(hpName,'(',leftright,')'));
        legend('2^n Order Harmonic Distortion','Other Harmonic Distortion','FontSize',15);
        ylabel('dBr');
        xlabel('Harmonic Orders');
        savePlot(path,'hd',LR);
        figure(figureNum);
    end
    
    subplot(5,2,10);
    
    if showFigure
        plot(thd3);
        ylim([0 0.003]);
        title(folderName);
    end
    
    fs3=fs2;
    freq3=log10(21):1/fs3:log10(1e4);
    thd3=mean(thd3);
    snr3=interp1(log10(hd(:,1)),hd(:,3),freq3);
    snr3=mean(db2mag(snr3));
    thdScore=thd3;
else
    thdScore=0;
    snrScore=0;
end
if figureNum==1
    result(ismember(result.name,{folderName}),:).THDL=thdScore;
else
    result(ismember(result.name,{folderName}),:).THDR=thdScore;
end

if isfile(mtdFile)
    mtdScore=csvread(mtdFile);
    mtdScore=mtdScore(1);
else
    mtdScore=0;
end
if figureNum==1
    result(ismember(result.name,{folderName}),:).MTDL=mtdScore;
else
    result(ismember(result.name,{folderName}),:).MTDR=mtdScore;
end

end

function y=phaseOffset(freq,phase,deltatime)
y=zeros(1,length(phase));
for i=1:length(phase)
    y(i)=phase(i)-deltatime*freq(1,i)*2*pi;
    y(i)=mod(y(i),2*pi);
    if y(i)>pi
        y(i)=y(i)-2*pi;
    end
end
end

function y=groupDelay(freq,phase)
phase=unwrap(phase);
y=ones(1,length(phase));
y(1)= -1*((phase(2)-phase(1))/(freq(2)-freq(1)))/(2*pi);
for n=2:1:length(phase)-1
    y(n)= -1*(((phase(n)-phase(n-1))/(freq(n)-freq(n-1)))+((phase(n+1)-phase(n))/(freq(n+1)-freq(n))))/(4*pi);
end
y(n)=-1*((phase(n)-phase(n-1))/(freq(n)-freq(n-1)))/(2*pi);
y=smooth(y,11);
end

function y=phaseDelay(freq,phase)
y=zeros(1,length(phase));
for n=1:length(phase)
    y(n)=phase(n)/freq(n)/(2*pi);
end
end

function csd=mycsd(y,fs,sliceinterval,windowinterval,starttime,endtime)
for k=fs*starttime:fs*sliceinterval:fs*endtime
    yy=y(k:min([k+(fs*windowinterval),length(y)]));
    ydft = fft(yy);
    ydft2 = ydft(1:(length(yy)-1)/2+1);
    freq = 0:fs/length(yy):fs/2;
    semilogx(freq,smooth(mag2db(abs(ydft2)),11));
    hold on;
end
hold off;
axis([20 20000 -90 0])
grid;
title('Decay');
xlabel('Hz');
ylabel('dB');
end

function [wavL,wavR,hdL,hdR,mtdL,mtdR,folderName,path]=readWav(path)
global hpName;
if isempty(path)
    path = uigetdir('..\5_2');
end
spl=strsplit(path,'\');
folderName=char(spl(end));
hpName=folderName;
path=join([path,'\']);
wavL=join([path,'l.wav']);
wavR=join([path,'r.wav']);
hdL=join([path,'\distortion\hdL.csv']);
hdR=join([path,'\distortion\hdR.csv']);
mtdL=join([path,'\distortion\mtdL.csv']);
mtdR=join([path,'\distortion\mtdR.csv']);
disp(path);
end

function phaseDiff=getPhaseDiff(phase2,phasem2,isFlipped,weightType)
global ptr10000hz;
if isFlipped
    k=pi;
else
    k=0;
end
phaseDiff=0;
for i=1:length(phase2)
    tmp=mod(abs(phase2(i)+k-phasem2(i)),2*pi);
    if tmp>pi
        tmp=abs(tmp-2*pi);
    end
    if weightType==1
        phaseDiff=phaseDiff+(2-2*i/length(phase2))*(tmp)/length(phase2);
    elseif weightType==2
        if i<=ptr10000hz
            m=2;
        else
            m=1;
        end
        phaseDiff=phaseDiff+m*(tmp)/(ptr10000hz*2+length(phase2)-ptr10000hz);
    end
end
end

function y=fs2ptr(fs2,myfreq)
y=round((log10(myfreq)-log10(20))*fs2);
end

