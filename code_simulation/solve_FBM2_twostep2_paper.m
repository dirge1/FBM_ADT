function L=solve_FBM2_twostep2_paper(pab,Data,Time,s,miu0)
%% ¶ÁÈ¡Êý¾Ý
r=0;

alpha_0=pab(1);

n=[];
for i=1:length(s)
    n=[n size(Data{i},2)];
end

for i=1:length(s)
    for j=1:n(i)
        a(i,j)=miu0(i,j)/exp(alpha_0*s(i));
    end
end
miu_a=mean(mean(a));
for i=1:length(s)
    for j=1:n(i)
        sa(i,j)=(miu0(i,j)/exp(alpha_0*s(i))-miu_a)^2;
    end
end
sigma_a=sqrt(mean(mean(sa)));

L=0;
for i=1:length(s)
    for j=1:n(i)
        L=L-1/2*(log(2*pi)+log(sigma_a^2*exp(2*alpha_0*s(i)))+(miu0(i,j)-miu_a*exp(alpha_0*s(i)))^2/(sigma_a^2*exp(2*alpha_0*s(i))));
    end
end

L=-L;
end