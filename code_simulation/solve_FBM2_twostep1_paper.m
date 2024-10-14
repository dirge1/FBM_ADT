function L=solve_FBM2_twostep1_paper(pab,Data,Time,s)
%% ¶ÁÈ¡Êý¾Ý

beita=pab(1);
H=pab(2);

n=[];
for i=1:length(s)
    n=[n size(Data{i},2)];
end

    P=0;
for i=1:length(s)
    t=Time{i}(2:end,1);
    T=t.^beita;
    N=size(Data{i},1)-1;
    P=P+N*size(Data{i},2);
    SIG=[];
    for p=1:N
        for q=1:N
            SIG(p,q)=(1/2*(t(p)^(2*H)+t(q)^(2*H)-(abs(t(p)-t(q)))^(2*H)));
        end
    end
    for j=1:n(i)
        Y0=Data{i}(2:end,j);
        miu0(i,j)=(Y0'*inv(SIG)*T)/(T'*inv(SIG)*T);
        sigma0(i,j)=(Y0-miu0(i,j)*T)'*inv(SIG)*(Y0-miu0(i,j)*T);
    end
end

sigma1=sqrt(sum(sum(sigma0))/P);

L=0;
for i=1:length(s)
    t=Time{i}(2:end,1);
    T=t.^beita;
    N=size(Data{i},1)-1;
    SIG=[];
    for p=1:N
        for q=1:N
            SIG(p,q)=1/2*(t(p)^(2*H)+t(q)^(2*H)-(abs(t(p)-t(q)))^(2*H));
        end
    end
    for j=1:n(i)
        Y=Data{i}(2:end,j);
        L=L-1/2*(N*log(2*pi)+log(det(SIG*sigma1^2))+(Y-miu0(i,j)*T)'*inv(SIG)*(Y-miu0(i,j)*T)/sigma1^2);
    end
end

L=-L;
end