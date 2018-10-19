Dreal = load('data_real.txt');
Dimag = load('data_imag.txt');

DataIQ = Dreal + 1j*Dimag ;

figure; imagesc(abs(DataIQ));


% vreal = load('ventana_real.txt');
% vimag = load('ventana_imag.txt');

% figure; plot(vreal);
% figure; plot(vimag);

te = 10e-6;
c = 3e8;
PRF = 1000;
fs = 50e6;
Tventana = 1e-5;
tau = 2.5e-6;
D = 1;
Fc = 5e9;
lambda = c/Fc;

vmax = PRF/2*lambda/2;

Rmin = te*c/2;
Rmax = Rmin + Tventana*c/2;

tita = -pi/2:pi/2000:pi/2;

E_tita = sin(pi*(D/lambda)*sin(tita))./(pi*(D/lambda)*sin(tita));

figure; plot(tita*180/pi, E_tita);
figure; plot(tita*180/pi, 10*log10(abs(E_tita).^2));




xdata = DataIQ(160,520:620);
v_index1 = -vmax:2*vmax/length(xdata):vmax-2*vmax/length(xdata);
figure; plot(v_index1,abs(fftshift(fft(xdata))));
xdata1 = DataIQ(400,520:680);
v_index2 = -vmax:2*vmax/length(xdata1):vmax-2*vmax/length(xdata1);
figure; plot(v_index2,abs(fftshift(fft(xdata1))));


