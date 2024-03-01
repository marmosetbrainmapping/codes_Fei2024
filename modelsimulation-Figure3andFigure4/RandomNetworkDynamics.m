function [V] = RandomNetworkDynamics(N,Nregion,W,Vrest,Vthres,r,beta,tau,dt,T,perturbations)


    V0 = unifrnd(Vrest,0,N*Nregion,1);
    V = zeros(N*Nregion,T);
    V(:,1) = V0;

    

    tic
    for t = 2:T

        Vprev = V(:,t-1);
        dV = 1/tau * ( r*Vprev + W*0.5*(1+erf(beta*(Vprev+Vthres))) + Vrest + perturbations(:,t) ) * dt;
        Vcurrent = Vprev + dV;
        V(:,t) = Vcurrent;

%         if mod(t,T/10) == 0
%             disp(['Complete-',num2str(t/T*100),'%'])
%             toc
%         end
    end


end