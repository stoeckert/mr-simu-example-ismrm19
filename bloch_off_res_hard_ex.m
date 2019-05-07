
% off-resonance exictation: 1 ms hard-pulse (b1 field on x-axis, no relaxation) 

F = linspace(-2*pi*4,2*pi*4,400);                           % offset frequencies [rad/ms] 
w = pi/2;                                                   % flipangle in radiant

%a) STA SOLUTION
M_STA = abs(sin(w)*sinc(F/(2*pi))); 

M_ODE=zeros(size(F)); M_ROT=zeros(size(F)); i=0;            % intialize result vectors
for freq = F; i=i+1;                                                      
%b) ODE SOLUTION
    bloch = @(t,M) [freq*M(2); -freq*M(1)+w*M(3); -w*M(2)];  
    [~,M] = ode45(bloch,[0 1],[0;0;1]);                     % integrate ODE
    M_ODE(i)=sqrt(M(end,1)^2+M(end,2)^2);                   % transv. magn.
%c) ROTATION SOLUTION 
    W = -sqrt(w^2+freq.^2);                                 % effective field rotation angle 
    n = [w 0 freq]/abs(W);                                  % effective field rotation axis
    a = cos(W/2)-1i*n(3)*sin(W/2);                          % Cayley-Klein Parameter (alpha)
    b = -1i*n(1)*sin(W/2);                                  % Cayley-Klein Parameter (beta)
    M_ROT(i)=abs(2*a'*b);                                   % transv. magn.
end

%plot result
F= F * 2*pi;
plot(F,M_STA,F(2:2:end),M_ODE(2:2:end),'ko',F(1:2:end),...
     M_ROT(1:2:end),'m','linewidth',2)
grid,legend({'STA','ODE','ROT'},'fontsize',14)
xlabel('frequency offset / kHz','fontsize',14)
ylabel('transverse magnetisation','fontsize',14)
title(sprintf('flip angle %2.0f degree',w*180/pi),'fontsize',15)
