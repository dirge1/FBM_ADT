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

            x0=[2 0.01 0.5 0.4 2e-6 4e-7];   %仿真初值

            alpha_0=x0(1);
            sigma=x0(2);
            beita=x0(3);
            H=0.5;
            miu_a=x0(5);
            sigma_a=x0(6);

            x00=[x0(3) x0(4)];
            lb = [0 0];
            ub = [3 1];

            A = [];
            b = [];
            Aeq = [];
            beq = [];
            fun=@(x) solve_FBM2_twostep1_paper(x,Data,Time,s);
            [pab, ~, flag1] = fmincon(fun,x00,A,b,Aeq,beq,lb,ub);

            % [pab, ~, flag1] = fmincon(@solve_FBM2_twostep1,x00,A,b,Aeq,beq,lb,ub);

            beita=pab(1);
            H=pab(2);

            P=0;

            for i=1:length(s0)
                t=Time0{i}(2:end,1);
                T=t.^beita;
                N=size(Data0{i},1)-1;
                P=P+N*size(Data0{i},2);
                SIG=[];
                for p=1:N
                    for q=1:N
                        SIG(p,q)=(1/2*(t(p)^(2*H)+t(q)^(2*H)-(abs(t(p)-t(q)))^(2*H)));
                    end
                end
                for j=1:n(i)
                    Y0=Data0{i}(2:end,j);
                    miu0(i,j)=(Y0'*inv(SIG)*T)/(T'*inv(SIG)*T);
                    sigma0(i,j)=(Y0-miu0(i,j)*T)'*inv(SIG)*(Y0-miu0(i,j)*T);
                end
            end

            sigma=sqrt(sum(sum(sigma0))/P);

            x00=[x0(1)];
            lb = [0];
            ub = [10];

            A = [];
            b = [];
            Aeq = [];
            beq = [];
            % [pab2, ~, flag1] = fmincon(@solve_FBM2_twostep2_paper,x00,A,b,Aeq,beq,lb,ub);
            fun=@(x) solve_FBM2_twostep2_paper(x,Data,Time,s,miu0);
            [pab2, ~, flag1] = fmincon(fun,x00,A,b,Aeq,beq,lb,ub);

            alpha_0=pab2(1);

            for i=1:length(s0)
                for j=1:n(i)
                    a(i,j)=miu0(i,j)/exp(alpha_0*s0(i));
                end
            end
            miu_a=mean(mean(a));
            for i=1:length(s0)
                for j=1:n(i)
                    sa(i,j)=(miu0(i,j)/exp(alpha_0*s0(i))-miu_a)^2;
                end
            end
            sigma_a=sqrt(mean(mean(sa)));

            x0=[alpha_0 sigma beita H miu_a sigma_a];


            alpha_0=x0(1);
            sigma=x0(2);
            beita=x0(3);
            H=0.5;
            miu_a=x0(5);
            sigma_a=x0(6);

            for o=1:10000
                for i=1:length(s0)
                    t=Time0{i}(2:end,1);
                    fai=exp(alpha_0*s0(i))*t.^beita;
                    N=size(Data0{i},1)-1;
                    SIG=[];
                    for p=1:N
                        for q=1:N
                            SIG(p,q)=(1/2*(t(p)^(2*H)+t(q)^(2*H)-(abs(t(p)-t(q)))^(2*H)));
                        end
                    end
                    for j=1:n(i)
                        Y0=Data0{i}(2:end,j);
                        miu0(i,j)=(Y0'*inv(SIG)*fai*sigma_a^2+miu_a*sigma^2)/(fai'*inv(SIG)*fai*sigma_a^2+sigma^2);
                        sigma0(i,j)=sqrt(sigma^2*sigma_a^2/(fai'*inv(SIG)*fai*sigma_a^2+sigma^2));
                    end
                end

                miu_a1=0;
                sigma_a1=0;

                for i=1:length(s0)
                    for j=1:n(i)
                        miu_a1=miu_a1+miu0(i,j);
                    end
                end
                miu_a1=miu_a1/sum(n);

                for i=1:length(s0)
                    for j=1:n(i)
                        sigma_a1=sigma_a1+miu0(i,j)^2+sigma0(i,j)^2-2*miu0(i,j)*miu_a1+miu_a1^2;
                    end
                end
                sigma_a1=sqrt(sigma_a1/sum(n));

                % x00=[x0(3) x0(4) x0(6)];
                x00=[x0(1) x0(3) x0(4)];
                lb = [0 0 0.5];
                ub = [10 3 0.5];

                A = [];
                b = [];
                Aeq = [];
                beq = [];

                % [pab, ~, flag1] = fmincon(@solve_FBM2_EM,x00,A,b,Aeq,beq,lb,ub);

                fun=@(x) solve_BM2_EM_paper(x,Data,Time,s,miu_a1,sigma_a1,miu0,sigma0);
                [pab, ~, flag1] = fmincon(fun,x00,A,b,Aeq,beq,lb,ub);

                alpha_01=pab(1);
                beita1=pab(2);
                H1=pab(3);

                sigma1=0;
                T=0;
                for i=1:length(s0)
                    t=Time0{i}(2:end,1);
                    fai=exp(alpha_01*s0(i))*t.^beita1;
                    N=size(Data0{i},1)-1;
                    SIG=[];
                    for p=1:N
                        for q=1:N
                            SIG(p,q)=(1/2*(t(p)^(2*H1)+t(q)^(2*H1)-(abs(t(p)-t(q)))^(2*H1)));
                        end
                    end
                    for j=1:n(i)
                        Y0=Data0{i}(2:end,j);
                        sigma1=sigma1+Y0'*inv(SIG)*Y0-2*miu0(i,j)*Y0'*inv(SIG)*fai+(miu0(i,j)^2+sigma0(i,j)^2)*fai'*inv(SIG)*fai;
                        T=T+length(t);
                    end
                end
                sigma1=sqrt(sigma1/T);

                T=[alpha_0 sigma beita H miu_a sigma_a];
                T1=[alpha_01 sigma1 beita1 H1 miu_a1 sigma_a1];
                epsilon(o)=sum(abs(T-T1)./T);
                alpha_0=alpha_01;
                beita=beita1;
                H=H1;
                sigma_a=sigma_a1;
                miu_a=miu_a1;
                sigma=sigma1;
                if  epsilon(o) < epsilon0
                    break;
                end
            end
            L=0;
            for i=1:length(s0)
                for j=1:n(i)
                    t=[];Y0=[];f=[];omiga=[];
                    t=Time0{i}(2:end,1);
                    Y0=Data0{i}(2:end,j);
                    f=miu_a*exp(alpha_0*s0(i))*t.^beita;
                    N=length(Y0);
                    for p=1:N
                        for q=1:N
                            omiga(p,q)=sigma^2*(1/2*(t(p)^(2*H)+t(q)^(2*H)-(abs(t(p)-t(q)))^(2*H)))+sigma_a^2*exp(2*alpha_0*s0(i))*t(p)^beita*t(q)^beita;
                        end
                    end
                    L=L-1/2*(N*log(2*pi)+log(det(omiga))+(Y0-f)'*inv(omiga)*(Y0-f));
                end
            end
            Q(ww)=L;
            AIC(ww)=-2*Q(ww)+2*5;
            result_two_step{pp,qq}(ww,:)=x0;
            result_EM{pp,qq}(ww,:)=T1;
        end
        %%
        pab_real=[2.5 0.1 1.5 0.1 1e-5 2e-6];
        mean_EM{pp,qq}=mean(result_EM{pp,qq});
        for i=1:6
            error_EM{pp,qq}(i)=abs(mean_EM{pp,qq}(i)-pab_real(i))/pab_real(i);
        end
        RE_EM{pp,qq}=sum(error_EM{pp,qq});

        mean_two_step{pp,qq}=mean(result_two_step{pp,qq});
        for i=1:6
            error_two_step{pp,qq}(i)=abs(mean_two_step{pp,qq}(i)-pab_real(i))/pab_real(i);
        end
        RE_two_step{pp,qq}=sum(error_two_step{pp,qq});

    end
end
