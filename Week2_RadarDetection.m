%% 

% Clean workspace
clear all; close all; clc


%% Fourier transform of sech(t)

L = 30;	% time slot to transform
n = 512;	% number of Fourier modes 2^9
t2 = linspace(-L,L,n+1);	% time discretization 
t = t2(1:n);		% only use the first n points (periodicity)
k = (2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; 	%frequency components (it looks
% weird because it has to correspond to the way DFTs are ordered in MATLAB)
% In the previous Week's code we reordered k from the get go, but here
% we'll need it for the filter later on.
u = sech(t);	%ideal signal in the time domain

% Plot the time domain
figure(1)
subplot(2,1,1)
plot(t,u, 'Linewidth', 2)
set(gca, 'Fontsize', 16)
xlabel('time (t)')
ylabel('|u|')

% Plot the frequency domain
ut = fft(u);
subplot(2,1,2)
plot(fftshift(k), fftshift(abs(ut)), 'r', 'Linewidth', 2)
set(gca, 'Fontsize', 16)
xlabel('frequency (k)')
ylabel('|ut|')

%% Adding noise

% Plot ideal signal
figure(2)
subplot(3,1,1)
plot(t,u,'Linewidth',2)
xlabel('time (t)')
ylabel('u')
set(gca,'Fontsize',16)

% Plot noisy signal (noise amplitude = 1)
noise = 1;
% Add white noise to the signal in frequency domain
utn = ut + noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));
un = ifft(utn); % Noisy signal in time domain
subplot(3,1,2)
plot(t,abs(un),'Linewidth',2)
xlabel('time (t)')
ylabel('u')
set(gca,'Fontsize',16)

% Plot noisy signal (noise amplitude = 5)
noise = 5;
% Add white noise to the signal in frequency domain
utn = ut + noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));
un = ifft(utn); % Noisy signal in time domain
subplot(3,1,3)
plot(t,abs(un),'Linewidth',2)
xlabel('time (t)')
ylabel('u')
set(gca,'Fontsize',16)

%% Even noisier signal

noise = 10;
utn = ut + noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));
un = ifft(utn);

% Plot the noisy signal 
figure(3)
subplot(2,1,1)
plot(t,abs(un),'Linewidth',2)
axis([-30 30 0 2])
xlabel('time (t)')
ylabel('|u|')
set(gca,'Fontsize',16)

% Plot the Fourier transform
subplot(2,1,2)
plot(fftshift(k), abs(fftshift(utn))/max(abs(fftshift(utn))),'r','Linewidth',2)
axis([-25 25 0 1])
xlabel('frequency (k)')
ylabel('|ut|/max(|ut|)')
set(gca,'Fontsize',16)

%% Frequency filtering with a Gaussian

tau = 0.2;
k0 = 0;
filter = exp(-tau*(k - k0).^2); % Define the filter
unft = filter.*utn; % Apply the filter to the signal in frequency space
unf = ifft(unft);

% Plot the unfiltered signal in the frequency domain and the Gaussian filter
figure(4)
subplot(3,1,1)
plot(fftshift(k),abs(fftshift(utn))/max(abs(fftshift(utn))),'r','Linewidth',2)
hold on
plot(fftshift(k),fftshift(filter),'k','Linewidth',2)
axis([-25 25 0 1])
xlabel('frequency (k)')
ylabel('|ut|/max(|ut|)')

% Plot the filtered signal in the frequency domain
subplot(3,1,2)
plot(fftshift(k),abs(fftshift(unft))/max(abs(fftshift(unft))),'r','Linewidth',2)
axis([-25 25 0 1])
xlabel('frequency (k)')
ylabel('|ut|/max(|ut|)')

% Plot the filtered signal in the time domain and the ideal signal
subplot(3,1,3)
plot(t,u,'k--','Linewidth',2)
hold on
plot(t,abs(unf),'b','Linewidth',2)
xlabel('time (t)')
ylabel('|u|')

%% No signal to filter

% What do we think will happen if there is nothing to filter?

utn = noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));
tau = 0.2;
k0 = 0;
filter = exp(-tau*(k - k0).^2); % Define the filter
unft = filter.*utn; % Apply the filter to the signal in frequency space
unf = ifft(unft);

% Plot the unfiltered signal in the frequency domain and the Gaussian
% filter
figure(5)
subplot(3,1,1)
plot(fftshift(k),abs(fftshift(utn))/max(abs(fftshift(utn))),'r','Linewidth',2)
hold on
plot(fftshift(k),fftshift(filter),'k','Linewidth',2)
axis([-25 25 0 1])
xlabel('frequency (k)')
ylabel('|ut|/max(|ut|)')

% Plot the filtered signal in the frequency domain
subplot(3,1,2)
plot(fftshift(k),abs(fftshift(unft))/max(abs(fftshift(unft))),'r','Linewidth',2)
axis([-25 25 0 1])
xlabel('frequency (k)')
ylabel('|ut|/max(|ut|)')

% Plot the filtered signal in the time domain and the ideal signal
subplot(3,1,3)
plot(t,u,'k--','Linewidth',2)
hold on
plot(t,abs(unf),'b','Linewidth',2)
xlabel('time (t)')
ylabel('|u|')

%% Frequency filtering around the wrong frequency

tau = 0.2;
k0 = 15;
utn = ut + noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));
filter = exp(-tau*(k - k0).^2); % Define the filter
unft = filter.*utn; % Apply the filter to the signal in frequency space
unf = ifft(unft);

% Plot the unfiltered signal in the frequency domain and the Gaussian
% filter
figure(6)
subplot(3,1,1)
plot(fftshift(k),abs(fftshift(utn))/max(abs(fftshift(utn))),'r','Linewidth',2)
hold on
plot(fftshift(k),fftshift(filter),'k','Linewidth',2)
axis([-25 25 0 1])
xlabel('frequency (k)')
ylabel('|ut|/max(|ut|)')

% Plot the filtered signal in the frequency domain
subplot(3,1,2)
plot(fftshift(k),abs(fftshift(unft))/max(abs(fftshift(unft))),'r','Linewidth',2)
axis([-25 25 0 1])
xlabel('frequency (k)')
ylabel('|ut|/max(|ut|)')

% Plot the filtered signal in the time domain and the ideal signal
subplot(3,1,3)
plot(t,u,'k--','Linewidth',2)
hold on
plot(t,abs(unf),'b','Linewidth',2)
xlabel('time (t)')
ylabel('|u|')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Creating a noisy signal

L = 30; % timeslot [-L,L]
n = 512; % number of Fourier modes

t2 = linspace(-L,L,n+1);
t = t2(1:n);
k = (2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; 
%define this so we don't need to use fftshift(k) every time we plot
ks = fftshift(k);
u = sech(t);
ut = fft(u);

% Add noise
noise = 10;
utn = ut + noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));

% Plot noisy signal
plot(ks,fftshift(abs(utn))/max(abs(utn)),'r','Linewidth',2)
set(gca,'Fontsize',16)
xlabel('frequency (k)')
ylabel('|ut|')

%% Average over 5 realizations of the signal

ave = zeros(1,n);
for j = 1:5
   utn = ut + noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));
   ave = ave + utn;

ave_plot = abs(fftshift(ave))/5;

% Plot averaged signal
plot(ks,ave_plot/max((ave_plot)),'r','Linewidth',2)
set(gca,'Fontsize',16)
xlabel('frequency (k)')
ylabel('|ut|')
pause
end

%% Comparing over different realizations

realize = [1 2 5 200];

for jj = 1:length(realize)
   u = sech(t);
   ave = zeros(1,n);
   ut = fft(u);
   
   % Averaging over realizations
   for j = 1:realize(jj)
      utn(j,:) = ut + noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));
      ave = ave+utn(j,:);
   end
   ave = abs(fftshift(ave))/realize(jj);
   
   % Plot averaged signals in one plot
   subplot(length(realize),1,jj)
   plot(ks,ave/max(ave),'r','Linewidth',2)
   set(gca,'Fontsize',16)
   axis([-20 20 0 1])
   ylabel('|fft(u)|','Fontsize',16)
end

hold on
plot(ks,abs(fftshift(ut))/max(abs(ut)),'k:','Linewidth',2)
xlabel('frequency (k)')

%% Shifted realizations

slice = 0:0.5:10;
[T,S] = meshgrid(t,slice);
[K,S] = meshgrid(k,slice);

U = sech(T - 10*sin(S)).*exp(1i*0*T);
figure(5)
subplot(2,1,1)
waterfall(T,2*S+1,U), colormap([0 0 0]), view(-15,70)
set(gca,'Fontsize',16,'Xlim',[-30 30],'Zlim',[0 2])
xlabel('time (t)'), ylabel('realizations'), zlabel('|u|')

for j = 1:length(slice)
   Ut(j,:) = fft(U(j,:));
   Kp(j,:) = fftshift(K(j,:));
   Utp(j,:) = fftshift(Ut(j,:));
   Utn(j,:) = Ut(j,:) + noise*(normrnd(0,1,1,n) + 1i*normrnd(0,1,1,n));
   Utnp(j,:) = fftshift(Utn(j,:))/max(abs(Utn(j,:)));
   Un(j,:) = ifft(Utn(j,:));
end

subplot(2,1,2)
waterfall(Kp,2*S+1,abs(Utp)/max(abs(Utp(1,:)))), colormap([0 0 0]), view(-15,70)
set(gca,'Fontsize',16,'Xlim',[-28 28])
xlabel('frequency (k)'), ylabel('realizations'), zlabel('|fft(u)|')

%% Noisy realizations

figure(6)
subplot(2,1,1)
waterfall(T,2*S+1,abs(Un)), colormap([0 0 0]), view(-15,70)
set(gca,'Fontsize',16,'Xlim',[-30 30],'Zlim',[0 2])
xlabel('time (t)'), ylabel('realizations'), zlabel('|u|')

subplot(2,1,2)
waterfall(Kp,2*S+1,abs(Utnp)/max(abs(Utnp(1,:)))), colormap([0 0 0]), view(-15,70)
set(gca,'Fontsize',16,'Xlim',[-28 28])
xlabel('time (t)'), ylabel('realizations'), zlabel('|fft(u)|')

%% Average over time and frequency domains to compare

Uave = zeros(1,n);
Utave = zeros(1,n);
for j = 1:length(slice)
   Uave = Uave + Un(j,:);
   Utave = Utave + Utn(j,:);
end
Uave = Uave/length(slice);
Utave = fftshift(Utave)/length(slice);

figure(7) 
subplot(2,1,1)
plot(t,abs(Uave),'b')
set(gca,'Fontsize',16)
xlabel('time (t)'), ylabel('|u|')

subplot(2,1,2)
plot(t,abs(Utave)/max(abs(Utave)),'r')
hold on
plot(ks,abs(fftshift(Ut(1,:))/max(abs(Ut(1,:)))),'k--','Linewidth',2)
axis([-20 20 0 1])
set(gca,'Fontsize',16)
xlabel('frequency (k)'), ylabel('|fft(u)|')
