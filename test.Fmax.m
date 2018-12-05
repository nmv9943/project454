
%Create Table
function [dp, dq, Pcom, Qcom] 

disp(' NewtonRaphson Table ');


disp('-----------------------------------------------------------------------------------------');
disp('| ybus |    Voltage  |  Theta  |     Injection    |   Generation   |       Load      |');
disp('| pu |           pu  |  Degree |     MW  |  MVAR  |   MW  |  MVAR  |  MW   | MVAR    | ');

for i=0
    Fmax i=20;
    
   
    fprintf('  %f', dp(i)); fprintf('   %f', dq(i)); %print 
    fprintf('  %f', pcom(i)); fprintf('   %f', Qcomp(i)); 
   
end

disp(' Line FLow Table');

disp('-------------------------------------------------------------------------------------');
disp('|From|To |    P    |    Q     | From| To |    P     |   Q     |      Line Loss      |');
disp('|Bus |Bus|   MW    |   MVar   | Bus | Bus|    MW    |  MVar   |     MW   |    MVar  |');

for i=0
    p=Q(i);
    
    
end





%NR method to calculate FMAX

i=0;
Fmax = 20;

epsilon = 1e-6;

estimationerror = 1;

while (estimationerror > epsilon && i< Fmax)
    
        i=i+1;%increment till error greater than epislon 
    
    S = V.p*conj(I);
    
    I=Y*V.p; V.p = V *1i*theta;
    
    Re=realpower;
    
    Ra=reactivepower;
    
    mismatch = [real(S(Re))-(p(Re)); imaginary(S(Ra))-(q(Ra))];
    


    


