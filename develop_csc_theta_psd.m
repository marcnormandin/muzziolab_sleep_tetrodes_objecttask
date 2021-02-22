figure
[Pxx,F] = pwelch(lowpass(cscz, 30,fs), [],[],[],fs);
plot(F, abs(Pxx))
grid on

%%
APxx = abs(Pxx);

fmin = 0;
fmax = 30;
f = linspace(fmin, fmax, 30);

x = [];
y = [];
for i = 1:(length(f)-1)
    ind = find(F >= f(i) & F <= f(i+1));
    s = mean(APxx(ind));
    x(i) = (f(i)+f(i+1))/2;
    y(i) = s;
end

figure
plot(x,y, 'r-')

%%
figure
plot(F, movmean(APxx, 100))
a = axis;
axis([0, 30, a(3), a(4)])
grid on
