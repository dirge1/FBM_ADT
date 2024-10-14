%% 考虑个体误差的FBM
rng(1107)
clear;clc; close all
pab_real=[2.5 0.1 1.5 0.1 1e-5 2e-6];
IN_index=30;  %检测次数
vec_no_index = 18; % 样本数
for pp=1:1
    for qq=1:1
        IN=IN_index(pp);  %检测次数
        vec_no=vec_no_index(qq); % 样本数

        for ww=1:1000
            ww
            Y=[];
            Data=[];
            % alpha_0 sigma beita H miu_a sigma_a
            pab=[2.5 0.1 1.5 0.1 1e-5 2e-6];
            s=(1/(40+273.15)-1./([80 100 120]+273.15))/(1/(40+273.15)-1/(120+273.15));
            H=pab(4);
            n = IN*100;

            for w=1:length(s)
                T = zeros(3,length(H),length(n)); % Time for Three Methods over different H's and n's

                cov_total = zeros(max(n)+1,max(n)+1,length(H));
                for f = 1:length(H)
                    c = (1/(max(n)))*[0.0001,1:(length(cov_total)-1)];
                    y = c.';
                    cov_total(:,:,f) = bsxfun(@(c,y) 0.5.*(c.^(2*H(f))+y.^(2*H(f))-abs(y-c).^(2*H(f))), c, y);
                end
                %%Computation for Different n and H
                for g = 1:length(n)
                    Z = normrnd(0,1,vec_no,n(g));

                    for f = 1:length(H)
                        cov = cov_total(1:n(g)+1,1:n(g)+1,f);

                        %%Cholesky Method
                        tic;
                        M = chol(cov,'lower');
                        B_chol{w} = M*[zeros(vec_no,1),Z]';
                        T(1,f,g) = toc;
                    end
                end
            end
            %%
            alpha_0=pab(1);
            miu_a=pab(5);
            sigma_a=pab(6);
            sigma=pab(2);
            beita=pab(3);
            a=normrnd(miu_a,sigma_a,[length(s),vec_no]);

            tt=(0:n)/n;
            T=n;
            for w=1:length(s)
                for i=1:vec_no
                    Y{w}(:,i)=a(w,i)*exp(alpha_0*s(w))*(tt*T).^beita+sigma*T^H.*B_chol{w}(:,i)';
                end
            end
            %% 选择数据
            r=0:100:n;
            for i=1:3
                for j=1:vec_no
                    Data{i}(:,j)=Y{i}(r+1,j);
                end
            end
            %%
            r=0:100:n;
            Time{1}=r';
            Time{2}=r';
            Time{3}=r';
            % save('simu_6_10.mat','Data','Time','H','alpha_0', 'sigma', 'beita', 'miu_a', 'sigma_a')

            r=0;
            epsilon0=1e-2;


            %% 参数估计
            q=0;
            for i=1:length(s)
                if i==r
                    q=1;
                else
                    Data0{i-q}=Data{i};
                    Time0{i-q}=Time{i};
                end
            end
            s0=s;
            if r>0
                s0(r)=[];
            end

            n=[];
            for i=1:length(s0)
                n=[n size(Data0{i},2)];
            end

            x0=[2 0.01 1 0.5 2e-6];   %仿真初值

            lb = [0 0 0 0.4 0];
            ub = [10 1 2 0.7 0.1];

            A = [];
            b = [];
            Aeq = [];
            beq = [];
            fun=@(x) solve_FBM_M0(x,Data,Time,s);
            options = optimoptions('fmincon', 'MaxFunctionEvaluations', 10000);
            [pab, ~, flag1(ww)] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub, [], options);

            % [pab, ~, flag1] = fmincon(@solve_FBM2_twostep1,x00,A,b,Aeq,beq,lb,ub);

            alpha_0=pab(1);

            sigma=pab(2);
            beita=pab(3);
            H=pab(4);
            miu_a=pab(5);
            L=0;
            for i=1:length(s0)
                for j=1:n(i)
                    t=[];Y=[];f=[];omiga=[];
                    t=Time0{i}(2:end,1);
                    Y=Data0{i}(2:end,j);
                    f=miu_a*exp(alpha_0*s(i))*t.^beita;
                    N=length(Y);
                    for p=1:N
                        for q=1:N
                            omiga(p,q)=sigma^2*(1/2*(t(p)^(2*H)+t(q)^(2*H)-(abs(t(p)-t(q)))^(2*H)));
                        end
                    end
                    L=L-1/2*(N*log(2*pi)+log(det(omiga))+(Y-f)'*inv(omiga)*(Y-f));
                end
            end
            Q(ww)=L;
            AIC(ww)=-2*Q(ww)+2*5;
            result_EM{pp,qq}(ww,:)=pab;
        end
        %%

        mean_EM{pp,qq}=mean(result_EM{pp,qq});

    end
end
