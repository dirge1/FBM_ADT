function L=solve_BM_M0(pab,Data,Time,s)
%% 读取数据
%电连接器应力松弛

%% 参数估计
n=[];
for i=1:length(s)
    n=[n size(Data{i},2)];
end
L=0;
alpha=pab(1);

sigma=pab(2);
beita=pab(3);
miu_a=pab(4);
H=0.5;
for i=1:length(s)
    for j=1:n(i)
        t=[];Y=[];f=[];omiga=[];
        t=Time{i}(2:end,1);
        Y=Data{i}(2:end,j);
        f=miu_a*exp(alpha*s(i))*t.^beita;
        N=length(Y);
        for p=1:N
            for q=1:N
                omiga(p,q)=sigma^2*(1/2*(t(p)^(2*H)+t(q)^(2*H)-(abs(t(p)-t(q)))^(2*H)));
            end
        end
        L=L-1/2*(N*log(2*pi)+log(det(omiga))+(Y-f)'*inv(omiga)*(Y-f));
    end
end
L=-L;
end