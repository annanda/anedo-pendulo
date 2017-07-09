' Serão feitos os modelos de acordo com Euler, RK de 2ª ordem e RK de 4ª ordem '

function z = gSh(Sh,Ih,Im)
    z = 40 - (0.00004 + (0.4/1000000)*Im)*Sh;
endfunction

function z = gIh(Sh,Ih,Im)
    z = ((0.8/1000000)*Im)*Sh - 0.1*Ih;
endfunction

function z = gIm(Ih,Im)
    z = ((0.4/1000000)*Ih)*(200000 - Im) - 0.25*Im;
endfunction

// Chm = 0.8 Bs = 0.4 ; 0.6 ; 0.8
// Cmh = 0.8 Bi = 0.8 ; 1.2 ; 1.6

// Resolução para Euler
function [t,Sh,Ih,Im] = Euler_Sistema(t0,tf,h,Sh0,Ih0,Im0)
    t = t0:h:tf
    n = length(t);
    Sh(1) = Sh0
    Ih(1) = Ih0
    Im(1) = Im0
    
    for i = 1:n-1
        KSh = gSh(Sh(i),Ih(i),Im(i))
        KIh = gIh(Sh(i),Ih(i),Im(i))
        KIm = gIm(Ih(i),Im(i))
        
        Sh(i+1) = Sh(i) + KSh*h;
        Ih(i+1) = Ih(i) + KIh*h;
        Im(i+1) = Im(i) + KIm*h;
    end
endfunction

// Resolução para RK2
function [t,Sh,Ih,Im] = RK2_Sistema(t0,tf,h,Sh0,Ih0,Im0)
    t = t0:h:tf
    n = length(t);
    Sh(1) = Sh0
    Ih(1) = Ih0
    Im(1) = Im0
    
    for i = 1:n-1
        K1Sh = gSh(Sh(i),Ih(i),Im(i))
        K1Ih = gIh(Sh(i),Ih(i),Im(i))
        K1Im = gIm(Ih(i),Im(i))
        
        K2Sh = gSh(Sh(i) + K1Sh*h/2,Ih(i) + K1Ih*h/2,Im(i) + K1Im*h/2)
        K2Ih = gIh(Sh(i) + K1Sh*h/2,Ih(i) + K1Ih*h/2,Im(i) + K1Im*h/2)
        K2Im = gIm(Ih(i) + K1Ih*h/2,Im(i) + K1Im*h/2)
        
        Sh(i+1) = Sh(i) + (K1Sh + K2Sh)*h/2;
        Ih(i+1) = Ih(i) + (K1Ih + K2Ih)*h/2;
        Im(i+1) = Im(i) + (K1Im + K2Im)*h/2;
    end
endfunction

function [t,Sh,Ih,Im] = RK23_Sistema(t0,tf,h,Sh0,Ih0,Im0)
    t = t0:h:tf
    n = length(t);
    Sh(1) = Sh0
    Ih(1) = Ih0
    Im(1) = Im0
    
    for i = 1:n-1
        K1Sh = gSh(Sh(i),Ih(i),Im(i))
        K1Ih = gIh(Sh(i),Ih(i),Im(i))
        K1Im = gIm(Ih(i),Im(i))
        
        K2Sh = gSh(Sh(i) + K1Sh*(3/4)*h,Ih(i) + K1Ih*(3/4)*h,Im(i) + K1Im*(3/4)*h)
        K2Ih = gIh(Sh(i) + K1Sh*(3/4)*h,Ih(i) + K1Ih*(3/4)*h,Im(i) + K1Im*(3/4)*h)
        K2Im = gIm(Ih(i) + K1Ih*(3/4)*h,Im(i) + K1Im*(3/4)*h)
        
        Sh(i+1) = Sh(i) + (K1Sh + 2*K2Sh)*h/3;
        Ih(i+1) = Ih(i) + (K1Ih + 2*K2Ih)*h/3;
        Im(i+1) = Im(i) + (K1Im + 2*K2Im)*h/3;
    end
endfunction

// Resolução para RK4
function [t,Sh,Ih,Im] = RK4_Sistema(t0,tf,h,Sh0,Ih0,Im0)
    t = t0:h:tf
    n = length(t);
    Sh(1) = Sh0
    Ih(1) = Ih0
    Im(1) = Im0
    
    for i = 1:n-1
        K1Sh = gSh(Sh(i),Ih(i),Im(i))
        K1Ih = gIh(Sh(i),Ih(i),Im(i))
        K1Im = gIm(Ih(i),Im(i))
        
        K2Sh = gSh(Sh(i) + K1Sh*h/2,Ih(i) + K1Ih*h/2,Im(i) + K1Im*h/2)
        K2Ih = gIh(Sh(i) + K1Sh*h/2,Ih(i) + K1Ih*h/2,Im(i) + K1Im*h/2)
        K2Im = gIm(Ih(i) + K1Ih*h/2,Im(i) + K1Im*h/2)
        
        K3Sh = gSh(Sh(i) + K2Sh*h/2,Ih(i) + K2Ih*h/2,Im(i) + K2Im*h/2)
        K3Ih = gIh(Sh(i) + K2Sh*h/2,Ih(i) + K2Ih*h/2,Im(i) + K2Im*h/2)
        K3Im = gIm(Ih(i) + K2Ih*h/2,Im(i) + K2Im*h/2)
        
        K4Sh = gSh(Sh(i) + K3Sh*h,Ih(i) + K3Ih*h,Im(i) + K3Im*h)
        K4Ih = gIh(Sh(i) + K3Sh*h,Ih(i) + K3Ih*h,Im(i) + K3Im*h)
        K4Im = gIm(Ih(i) + K3Ih*h,Im(i) + K3Im*h)
        
        Sh(i+1) = Sh(i) + (K1Sh + 2*K2Sh + 2*K3Sh + K4Sh)*h/6;
        Ih(i+1) = Ih(i) + (K1Ih + 2*K2Ih + 2*K3Ih + K4Ih)*h/6;
        Im(i+1) = Im(i) + (K1Im + 2*K2Im + 2*K3Im + K4Im)*h/6;
    end
endfunction

// Resolução para Dormand-Prince
function [t,Sh,Ih,Im] = ODE45_Sistema(t0,tf,h0,Sh0,Ih0,Im0)
    Sh(1) = Sh0
    Ih(1) = Ih0
    Im(1) = Im0
    t(1) = t0
    h(1) = h0
    
    while t($) < tf
        K1Sh = h($)*gSh(Sh($),Ih($),Im($))
        K1Ih = h($)*gIh(Sh($),Ih($),Im($))
        K1Im = h($)*gIm(Ih($),Im($))
        
        K2Sh = h($)*gSh(Sh($) + K1Sh/5,Ih($) + K1Ih/5,Im($) + K1Im/5)
        K2Ih = h($)*gIh(Sh($) + K1Sh/5,Ih($) + K1Ih/5,Im($) + K1Im/5)
        K2Im = h($)*gIm(Ih($) + K1Ih/5,Im($) + K1Im/5)
        
        K3Sh = h($)*gSh(Sh($) + K1Sh*(3/40) + K2Sh*(9/40),Ih($) + K1Ih*(3/40) + K2Ih*(9/40),Im($) + K1Im*(3/40) + K2Im*(9/40))
        K3Ih = h($)*gIh(Sh($) + K1Sh*(3/40) + K2Sh*(9/40),Ih($) + K1Ih*(3/40) + K2Ih*(9/40),Im($) + K1Im*(3/40) + K2Im*(9/40))
        K3Im = h($)*gIm(Ih($) + K1Ih*(3/40) + K2Ih*(9/40),Im($) + K1Im*(3/40) + K2Im*(9/40))
        
        K4Sh = h($)*gSh(Sh($) + K1Sh*(44/45) - K2Sh*(56/15) + K3Sh*(32/9),Ih($) + K1Ih*(44/45) - K2Ih*(56/15) + K3Ih*(32/9),Im($) + K1Im*(44/45) - K2Im*(56/15) + K3Im*(32/9))
        K4Ih = h($)*gIh(Sh($) + K1Sh*(44/45) - K2Sh*(56/15) + K3Sh*(32/9),Ih($) + K1Ih*(44/45) - K2Ih*(56/15) + K3Ih*(32/9),Im($) + K1Im*(44/45) - K2Im*(56/15) + K3Im*(32/9))
        K4Im = h($)*gIm(Ih($) + K1Ih*(44/45) - K2Ih*(56/15) + K3Ih*(32/9),Im($) + K1Im*(44/45) - K2Im*(56/15) + K3Im*(32/9))
        
        K5Sh = h($)*gSh(Sh($) + K1Sh*(19372/6561) - K2Sh*(25360/2187) + K3Sh*(64448/6561) - K4Sh*(212/729),Ih($) + K1Ih*(19372/6561) - K2Ih*(25360/2187) + K3Ih*(64448/6561) - K4Ih*(212/729),Im($) + K1Im*(19372/6561) - K2Im*(25360/2187) + K3Im*(64448/6561) - K4Im*(212/729))
        K5Ih = h($)*gIh(Sh($) + K1Sh*(19372/6561) - K2Sh*(25360/2187) + K3Sh*(64448/6561) - K4Sh*(212/729),Ih($) + K1Ih*(19372/6561) - K2Ih*(25360/2187) + K3Ih*(64448/6561) - K4Ih*(212/729),Im($) + K1Im*(19372/6561) - K2Im*(25360/2187) + K3Im*(64448/6561) - K4Im*(212/729))
        K5Im = h($)*gIm(Ih($) + K1Ih*(19372/6561) - K2Ih*(25360/2187) + K3Ih*(64448/6561) - K4Ih*(212/729),Im($) + K1Im*(19372/6561) - K2Im*(25360/2187) + K3Im*(64448/6561) - K4Im*(212/729))
        
        K6Sh = h($)*gSh(Sh($) + K1Sh*(9017/3168) - K2Sh*(355/33) - K3Sh*(46732/5247) + K4Sh*(49/176) - K5Sh*(5103/18656),Ih($) + K1Ih*(9017/3168) - K2Ih*(355/33) - K3Ih*(46732/5247) + K4Ih*(49/176) - K5Ih*(5103/18656),Im($) + K1Im*(9017/3168) - K2Im*(355/33) - K3Im*(46732/5247) + K4Im*(49/176) - K5Im*(5103/18656))
        K6Ih = h($)*gIh(Sh($) + K1Sh*(9017/3168) - K2Sh*(355/33) - K3Sh*(46732/5247) + K4Sh*(49/176) - K5Sh*(5103/18656),Ih($) + K1Ih*(9017/3168) - K2Ih*(355/33) - K3Ih*(46732/5247) + K4Ih*(49/176) - K5Ih*(5103/18656),Im($) + K1Im*(9017/3168) - K2Im*(355/33) - K3Im*(46732/5247) + K4Im*(49/176) - K5Im*(5103/18656))
        K6Im = h($)*gIm(Ih($) + K1Ih*(9017/3168) - K2Ih*(355/33) - K3Ih*(46732/5247) + K4Ih*(49/176) - K5Ih*(5103/18656),Im($) + K1Im*(9017/3168) - K2Im*(355/33) - K3Im*(46732/5247) + K4Im*(49/176) - K5Im*(5103/18656))
        
        K7Sh = h($)*gSh(Sh($) + K1Sh*(35/384) + K3Sh*(500/1113) + K4Sh*(125/192) - K5Sh*(2187/6784) + K6Sh*(11/84),Ih($) + K1Ih*(35/384) + K3Ih*(500/1113) + K4Ih*(125/192) - K5Ih*(2187/6784) + K6Ih*(11/84),Im($) + K1Im*(35/384) + K3Im*(500/1113) + K4Im*(125/192) - K5Im*(2187/6784) + K6Im*(11/84))
        K7Ih = h($)*gIh(Sh($) + K1Sh*(35/384) + K3Sh*(500/1113) + K4Sh*(125/192) - K5Sh*(2187/6784) + K6Sh*(11/84),Ih($) + K1Ih*(35/384) + K3Ih*(500/1113) + K4Ih*(125/192) - K5Ih*(2187/6784) + K6Ih*(11/84),Im($) + K1Im*(35/384) + K3Im*(500/1113) + K4Im*(125/192) - K5Im*(2187/6784) + K6Im*(11/84))
        K7Im = h($)*gIm(Ih($) + K1Ih*(35/384) + K3Ih*(500/1113) + K4Ih*(125/192) - K5Ih*(2187/6784) + K6Ih*(11/84),Im($) + K1Im*(35/384) + K3Im*(500/1113) + K4Im*(125/192) - K5Im*(2187/6784) + K6Im*(11/84))
        
        
        YSh = Sh($) + (35/384)*K1Sh + (500/1113)*K3Sh + (125/192)*K4Sh - (2187/6784)*K5Sh + (11/84)*K6Sh;
        YIh = Ih($) + (35/384)*K1Ih + (500/1113)*K3Ih + (125/192)*K4Ih - (2187/6784)*K5Ih + (11/84)*K6Ih;
        YIm = Im($) + (35/384)*K1Im + (500/1113)*K3Im + (125/192)*K4Im - (2187/6784)*K5Im + (11/84)*K6Im;
        
        ZSh = Sh($) + (5179/57600)*K1Sh + (7571/16695)*K3Sh + (393/640)*K4Sh - (92097/339200)*K5Sh + (187/2100)*K6Sh + (1/40)*K7Sh;
        ZIh = Ih($) + (5179/57600)*K1Ih + (7571/16695)*K3Ih + (393/640)*K4Ih - (92097/339200)*K5Ih + (187/2100)*K6Ih + (1/40)*K7Ih;
        ZIm = Im($) + (5179/57600)*K1Im + (7571/16695)*K3Im + (393/640)*K4Im - (92097/339200)*K5Im + (187/2100)*K6Im + (1/40)*K7Im;
        
        
        sSh = (h($)*(10^-5)/(2*(tf-t0)*abs(YSh - ZSh)))^(1/4);
        sIh = (h($)*(10^-5)/(2*(tf-t0)*abs(YIh - ZIh)))^(1/4);
        sIm = (h($)*(10^-5)/(2*(tf-t0)*abs(YIm - ZIm)))^(1/4);
        
        if min(sSh,sIh,sIm) >= 2
            if (tf - t($)) > 2*h($)
                h($+1) = 2*h($);
                t($+1) = t($) + h($);
            end
            if (tf - t($)) <= 2*h($)
                h($+1) = (tf - t($));
                t($+1) = tf;
            end 
            Sh($+1) = YSh;
        end
        
        if min(sSh,sIh,sIm) >= 1 & min(sSh,sIh,sIm) < 2
            if (tf - t($)) > h($)
                h($+1) = h($);
                t($+1) = t($) + h($);
            end
            if (tf - t($)) <= h($)
                h($+1) = (tf - t($));
                t($+1) = tf;
            end
            Sh($+1) = YSh;
        end
        
        
        if min([sSh,sIh,sIm]) < 1
            if (tf - t($)) > h($)/2
                h($) = h($)/2;
                t($) = t($) - h($);
            end
            if (tf - t($)) <= h($)/2
                h($) = (tf - t($));
                t($) = t($-1) + h($);
            end
        end
        
    end
    
endfunction

//[t1,Sh1,Ih1,Im1] = Euler_Sistema(0,365,0.1,1000000,1,0)
//[t2,Sh2,Ih2,Im2] = RK2_Sistema(0,365,0.1,1000000,1,0)
[t3,Sh3,Ih3,Im3] = RK4_Sistema(0,365,0.1,1000000,10,0)
//[t4,Sh4,Ih4,Im4] = RK23_Sistema(0,365,0.1,1000000,1,0)



f3 = figure()
f3.background = 8
f3.figure_name = 'Resolução por Runge-Kutta de 4ª Ordem com h = 0.1'
plot(t3,Sh3,'r:')
plot(t3,Ih3,'r-')
plot(t2,Im3,'r--')

