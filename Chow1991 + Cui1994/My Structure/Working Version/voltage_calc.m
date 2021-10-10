function V = voltage_calc(xxx,zzz,mode)
lex= 1*1e-3; lem= 1*1e-3; gap= 0; h=0.17*1e-3; te = 0.06*1e-3;
N1 = 10; N2 = 6; N_MOM= N1*N2;
dx=0.1*1e-3; dz=0.01*1e-3; %mm
eps0=8.85*1e-12;
epsp = 1; epsr = 2.7;
%% -- loading data
noPlasma = 0;
withPlasma = 1;

if strcmpi (mode, 'noPlasma')
    chosen_mode = noPlasma;
elseif strcmpi (mode, 'withPlasma')
    chosen_mode = withPlasma;
else
    disp ('ERROR: error in parsing the argument "mode".');
    return;
end
switch chosen_mode
    case noPlasma
        %epsp = 1;
        a1 = struct2array(load('a1_origGF_V03_h_edited_3.mat'));
        %a2 = struct2array(load('a1_new.mat'));
    case withPlasma
        %epsp = 1+77.62*1j;
        a1 = struct2array(load('a1_new.mat'));
        %a2 = struct2array(load('a1_new.mat'));
end
%% -- finding voltage
if  (((0 <= xxx) && (xxx <=lex)) && ((h <= zzz) && (zzz <= (h+te))))
   V = 6*1e3; 
else
%V_vect = zeros(2*N_MOM,1);
V_vect = zeros(N_MOM,1);
    m1 = 0;
    for xm1= 0.05*1e-3 : dx: (lex - 0.05*1e-3)
        for zm1= (h+0.005*1e-3): dz: (h+te - 0.005*1e-3)
            m1 = m1 + 1;
            G = GreenF_up_electrode(xxx,zzz,xm1,zm1, dx, dz, epsp, epsr);
            V_vect(m1,1) = real((a1(m1,1)) * (G/1)) ;
        end   
    end
% disp(V_vect)    
% m2 = 0;
%     for xm2= (lex + gap + 0.05*1e-3) : dx: (lex+ gap + lem - 0.05*1e-3)
%         for zm2= (0.005*1e-3): dz: (te - 0.005*1e-3)
%             m2 = m2 + 1;
%             V_vect(N_MOM+m2,1) = a2(m2,1) * GreenF(x,z,xm2,zm2, dx, dz);
%         end   
%     end
    V = sum(V_vect);
end
end