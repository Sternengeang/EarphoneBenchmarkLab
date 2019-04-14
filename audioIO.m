function yk=audioIO(xk,fs)
recObj = audiorecorder(fs,24,2);
record(recObj);
player=audioplayer(xk,fs,24);
play(player);
pause(numel(xk)/fs+1);
stop(recObj);
yk = getaudiodata(recObj);
yk=yk(:,2);
end