function [ydft, freq]=mydft(yk,fs)
ydft = fft(hanning(length(yk)).*yk);
ydft= ydft(1:(length(yk))/2+1);
freq = 0:fs/length(yk):fs/2;
end