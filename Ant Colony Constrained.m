%%
clc;
clear;
maximum_iterations = 1000;
number_of_variables = 5;
number_of_division = 100;
rho = 0.5;
eta = 2;
alpha = 1;
lower_bound = [78;33;27;27;27];
upper_bound =  [102;45;45;45;45];
number_of_ants = 276;
ants = zeros(number_of_variables,number_of_ants);
x = zeros(number_of_variables,number_of_division+1);
    for i = 1:number_of_division+1
        x(:,i) = lower_bound + (i-1)*(upper_bound-lower_bound)/number_of_division;
    end
tprob = zeros(number_of_variables,number_of_division+1);
tau = ones(number_of_variables,number_of_division+1);
newtau = ones(number_of_variables,number_of_division+1);
tx = zeros(number_of_variables,number_of_division+1);
itr = 0;
flag2 = 1;
A = [];
K = 1;
%% Code for simple penalty functioon based on l2 norm
    while (K <= 10E9)
        while(flag2)

        flag = 1;
        tausum = sum(tau,2);
        new_tau = tau.^alpha;
        new_tausum = sum(new_tau,2);
        for i = 1:number_of_variables
        tprob(i,:) = new_tau(i,:)/new_tausum(i);
        end

        tx = tprob;
        for i = 2:number_of_division+1
            tx(:,i) = tx(:,i) + tx(:,i-1);
        end

        rnd = rand(number_of_variables,number_of_ants);

        for k = 1:number_of_variables
            for j = 1:number_of_ants
                pos = find((minus(tx(k,:) , rnd(k,j)))>0 ,1);
                if (size(pos) == [1 0])
                    pos = number_of_division+1; %Last element
                end
            ants(k,j) = x(k,pos);
            end
        end

    %f = ants(1,:).^2 - 2.*ants(1,:) - 11;
    %RosenBrock function
    %f = (1-ants(1,:)).^2 + 100.*(ants(2,:) - ants(1,:).^2).^2;
    %f = ants(1,:).^2 + ants(2,:).^2;

    % Constraints
    c = struct('constraints',{});
    c(1).constraints = -0.0022053*ants(3,:).*ants(5,:) + 0.0056858*ants(2,:).*ants(5,:) + 0.0006262*ants(1,:).*ants(4,:) - 6.665593;
    c(2).constraints = 0.0022053*ants(3,:).*ants(5,:) - 0.0056858*ants(2,:).*ants(5,:) - 0.0006262*ants(1,:).*ants(4,:) - 85.334407;
    c(3).constraints = 0.0071317*ants(2,:).*ants(5,:) + 0.0021813*ants(3,:).*ants(3,:) + 0.0029955*ants(1,:).*ants(2,:) - 29.48751;
    c(4).constraints = -0.0071317*ants(2,:).*ants(5,:) - 0.0021813*ants(3,:).*ants(3,:) - 0.0029955*ants(1,:).*ants(2,:) + 9.48751;
    c(5).constraints = 0.0047026*ants(3,:).*ants(5,:) + 0.0019085*ants(3,:).*ants(4,:) + 0.0012547*ants(1,:).*ants(3,:) - 15.699039;
    c(6).constraints = -0.004702*ants(3,:).*ants(5,:) - 0.0019085*ants(3,:).*ants(4,:) - 0.0012547*ants(1,:).*ants(3,:) + 10.699039;

    f = [];res = zeros(1,number_of_ants);
    for i = 1:number_of_ants
        res = 0;
        for j = 1:size(c,2)
            penalty(j) = min(0,c(j).constraints(i));
        end
        res(1,i) = sum(penalty);
    end

    % Objective Function
    f = 37.293239*ants(1,:) + 0.8356891*ants(1,:).*ants(5,:)  + 5.3578547*ants(3,:).*ants(3,:) - 40792.141 + (K*res(1,:));

    best1 = [];best2 = [];best3 = [];best4 = [];best5 = [];
    [fbest,I] = min(f);                                                        
    best1 = [best1 find( x(1,:) == ants(1,I))];                               
    best2 = [best2 find( x(2,:) == ants(2,I))];
    best3 = [best3 find( x(3,:) == ants(3,I))];
    best4 = [best4 find( x(4,:) == ants(4,I))];
    best5 = [best5 find( x(5,:) == ants(5,I))];
    best = [best1;best2;best3;best4;best5];
    fworst = max(f);

        for i = 1:number_of_variables
        newtau(i,:) = (1-rho).*tau(i,:);
            for j =1:length(best(i,:))
                newtau(i,best(i,j)) = tau(i,best(i,j)) + length(best(i,:))*eta*fbest/fworst;
            end
            tau(i,:) = newtau(i,:);
        end

        for i = 1:number_of_ants
        if(ants(:,1) == ants(:,i))
            flag = flag + 1;
        end
        end

        if (flag == number_of_ants+1)
            flag2 = 0;
        end

        itr = itr + 1;

        for i=1:number_of_ants
            plot(ants(1,i),ants(2,i),'r+');
            hold on;
        end
        hold off;
        pause

        end
    flag2 = 1;
    K = K*10;
    end
    