function [deltatime, isFlipped]=getPeakDelta(y,fs)
[maxy, maxypos]=max(abs(y));
isFlipped=~(maxy==y(maxypos));
deltatime=-(maxypos-1)/fs;
end