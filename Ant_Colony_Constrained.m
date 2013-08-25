%%
clc;
clear;
maximum_iterations = 1000;
number_of_variables = 13;
number_of_division = 100;
rho = 0.5;
eta = 2;
alpha = 1;
lower_bound = [0;0;0;0;0;0;0;0;0;0;0;0;0];
upper_bound = [1;1;1;1;1;1;1;1;1;10E9;10E9;10E9;1];
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
%% Code for simple penalty function based on l2 norm
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
    c(1).constraints = 2*ants(1,:) + 2*ants(2,:) + ants(10,:) + ants(11,:) - 10;
    c(2).constraints = 2*ants(1,:) + 2*ants(3,:) + ants(10,:) + ants(12,:) - 10;
    c(3).constraints = 2*ants(2,:) + 2*ants(3,:) + ants(11,:) + ants(12,:) - 10;
    c(4).constraints = -8*ants(1,:) + ants(10,:);
    c(5).constraints = -8*ants(2,:) + ants(11,:);
    c(6).constraints = -8*ants(3,:) + ants(12,:);
    c(7).constraints = -2*ants(4,:) - ants(5,:) + ants(10,:);
    c(8).constraints = -2*ants(6,:) - ants(7,:) + ants(11,:);
    c(9).constraints = -2*ants(8,:) - ants(9,:) + ants(12,:);

    f = [];res = zeros(1,number_of_ants);
    for i = 1:number_of_ants
        res = 0;
        for j = 1:size(c,2)
            penalty(j) = min(0,c(j).constraints(i));
        end
        res(1,i) = sum(penalty);
    end

    % Objective Function
    f = 5*ants(1,:) + 5*ants(2,:) + 5*ants(3,:) + 5*ants(4,:) ...
        - 5*ants(1,:).*ants(1,:) - 5*ants(2,:).*ants(2,:) ...
        - 5*ants(3,:).*ants(3,:) - 5*ants(4,:).*ants(4,:) ...
        - ants(5,:) - ants(6,:) - ants(7,:) - ants(8,:) - ants(9,:) ...
        - ants(10,:) - ants(11,:) - ants(12,:) - ants(13,:) + (K*res(1,:));

    best1 = [];best2 = [];best3 = [];best4 = [];best5 = [];best6 = [];best7 = [];best8 = [];best9 = [];best10 = [];best11 = [];best12 = [];best13 = [];
    [fbest,I] = min(f);                                                        
    best1 = find( x(1,:) == ants(1,I));
    best2 = find( x(2,:) == ants(2,I));
    best3 = find( x(3,:) == ants(3,I));
    best4 = find( x(4,:) == ants(4,I));
    best5 = find( x(5,:) == ants(5,I));
    best6 = find( x(6,:) == ants(6,I));
    best7 = find( x(7,:) == ants(7,I));
    best8 = find( x(8,:) == ants(8,I));
    best9 = find( x(9,:) == ants(9,I));
    best10 = find( x(10,:) == ants(10,I));
    best11 = find( x(11,:) == ants(11,I));
    best12 = find( x(12,:) == ants(12,I));
    best13 = find( x(13,:) == ants(13,I));
    best = [best1;best2;best3;best4;best5;best6;best7;best8;best9;best10;best11;best12;best13];
    fworst = max(f);

        for i = 1:number_of_variables
        newtau(i,:) = (1-rho).*tau(i,:);
            for j =1:length(best(i,:))
                newtau(i,best(i,j)) = tau(i,best(i,j)) + (length(best(i,:))*eta*fbest/fworst);
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
            plot(ants(1,i),ants(2,i),'ro');
            hold on;
        end
        hold off;
        pause

        end
    flag2 = 1;
    K = K*10;
    end
    