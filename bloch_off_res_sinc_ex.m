
T     = linspace(0,5,200);                                                 % timesteps axis [ms]
dt    = T(2)-T(1);                                                         % sampling interval
F     = linspace(-2*pi*4,2*pi*4,200);                                      % offset frequencies [rad/ms] 
w     = pi/2;                                                              % tip-angle [rad]
shape = @(t) sinc(2*(t-2.5));                                              % sinc RF pulse 
v     = w / abs(trapz(T,shape(T)));                                        % normalized tip-angle

%a) STA SOLUTION using fft
N     = round(numel(F)*.5*pi/dt/max(F))*2;                                 % nmuber of points on fft
M_STA = sin(w)*abs(fftshift(fft(shape(T),N))) * dt*v/w; 
M_STA = M_STA(N/2+(-numel(F)/2+1:numel(F)/2));                             % extract frequencies of interest

M_ODE=zeros(size(F)); M_ROT=zeros(size(F)); i=0;                           % intialize result vectors
for freq = F; i=i+1;                                                      
%b) ODE SOLUTION                                      
  bloch = @(t,M) [ freq*M(2);-freq*M(1)+v*shape(t)*M(3);-v*shape(t)*M(2)];
  [~,M] = ode45(bloch,T,[0;0;1]);                                          % integrate ODE  
  M_ODE(i)=norm(M(end,1:2));                                               % transv. magn.
%c) ROTATION MATRIX SOLUTION 
  Q=eye(2);                                                                % init Cayley-Klein matrix
  for t=T
        W = -dt*sqrt(abs(v*shape(t))^2+freq^2);                            % effective field rotation angle
        n = dt * [v*shape(t) 0 freq]/abs(W);                               % effective field rotation axis
        a = cos(W/2)-1i*n(3)*sin(W/2);                                     % Cayley-Klein Parameter (alpha)
        b = (-1i*n(1)+n(2))*sin(W/2);                                      % Cayley-Klein Parameter (beta)
        Q=[a -b';b a']*Q;                                                  % chained rotation 
  end
  M_ROT(i) = abs(2*Q(1,1)'*Q(2,1));                                        % transv. magn.
end

%plot result
F=F/(2*pi);
plot(F,M_STA,F,M_ROT,'m',F(1:2:end),M_ODE(1:2:end),'ko'...
     ,'linewidth',2),axis([-4 4 0 1.2])
grid,legend({'STA','ROT','ODE'},'fontsize',14)
xlabel('frequency offset / kHz','fontsize',14)
ylabel('transverse magnetisation','fontsize',14)
title(sprintf('flip angle %2.0f degree',w*180/pi),'fontsize',15)


