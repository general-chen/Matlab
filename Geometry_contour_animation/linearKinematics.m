function [s] = linearKinematics(Uo,Ts,a,sNorm,c,t)
 
N = round(t/Ts);
sVec = zeros(N,1);
s = zeros(N,1);
i = 1;
sTot = 0;
 
while i <= N;
    if i == 1
    sVec(i) = Uo*Ts + 0.5*a*Ts^2;
    U = Uo + a*Ts;
    else
    sVec(i) = U*Ts + 0.5*a*Ts^2;
    U = U + a*Ts;
    end
     
    sTot = sTot +sVec(i)/c;
     
    if i == N && sTot/c ~= sNorm
        fprintf('Predicted distance = %0.2f \n',sTot)
    end
    s(i) = sTot;
    i = i+1;
end