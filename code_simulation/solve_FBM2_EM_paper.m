function L=solve_FBM2_EM_paper(pab,Data,Time,s,miu_a1,sigma_a1,miu0,sigma0)
%% ¶ÁÈ¡Êý¾Ý

alpha_0=pab(1);
beita=pab(2);
H=pab(3);

n=[];
for i=1:length(s)
    n=[n size(Data{i},2)];
end

sigma1=0;
T=0;
for i=1:length(s)
    t=Time{i}(2:end,1);
    fai=exp(alpha_0*s(i))*t.^beita;
    N=size(Data{i},1)-1;
    SIG=[];
    for p=1:N
        for q=1:N
            SIG(p,q)=(1/2*(t(p)^(2*H)+t(q)^(2*H)-(abs(t(p)-t(q)))^(2*H)));
        end
    end
    for j=1:n(i)
        Y0=Data{i}(2:end,j);
        % equation 31
        sigma1=sigma1+Y0'*inv(SIG)*Y0-2*miu0(i,j)*Y0'*inv(SIG)*fai+(miu0(i,j)^2+sigma0(i,j)^2)*fai'*inv(SIG)*fai;
        T=T+length(t);
    end
end
sigma1=sqrt(sigma1/T);

L=0;
for i=1:length(s)
    t=Time{i}(2:end,1);
    fai=exp(alpha_0*s(i))*t.^beita;
    N=size(Data{i},1)-1;
    Q=[];
    for p=1:N
        for q=1:N
            Q(p,q)=sigma1^2*(1/2*(t(p)^(2*H)+t(q)^(2*H)-(abs(t(p)-t(q)))^(2*H)));
        end
    end
    SIG=Q/sigma1^2;
    for j=1:n(i)
        Y=Data{i}(2:end,j);
        % equation 28
        L=L-1/2*((N+1)*log(2*pi)...
            +log(det(SIG))...
            +log(sigma_a1^2)...
            +(miu0(i,j)^2+sigma0(i,j)^2-2*miu0(i,j)*miu_a1+miu_a1^2)/sigma_a1^2+...
            +N*log(sigma1^2)...
            +(Y'*inv(SIG)*Y-2*miu0(i,j)*Y'*inv(SIG)*fai+(miu0(i,j)^2+sigma0(i,j)^2)*fai'*inv(SIG)*fai)/sigma1^2);
    end
end

L=-L;
end