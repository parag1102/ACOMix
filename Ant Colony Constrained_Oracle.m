%%
clc;
clear;
maximum_iterations = 1000;
number_of_variables = 5;
number_of_division = 100;
rho = 0.5;
eta = 2;
alpha2 = 2;
lower_bound = [78;33;27;27;27];                                             %Upper bound and Lower bound may be modified according to the 
upper_bound =  [102;45;45;45;45];                                           % constraints of the problem
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
omega = 10E8;                                                               % Since the Global Optimal Value is unknown, omega is assumed to be a large number
itr2 = 0;
alpha1 = 0;
%%
    while(flag2)                                                            %Same termination criterion for the oracle
    while(flag2)

    flag = 1;
    tausum = sum(tau,2);
    new_tau = tau.^alpha2;
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
    c = struct('constraints',{});                                           %Constraints written in the form of a structure
    c(1).constraints = -0.0022053*ants(3,:).*ants(5,:) + 0.0056858*ants(2,:).*ants(5,:) + 0.0006262*ants(1,:).*ants(4,:) - 6.665593;
    c(2).constraints = 0.0022053*ants(3,:).*ants(5,:) - 0.0056858*ants(2,:).*ants(5,:) - 0.0006262*ants(1,:).*ants(4,:) - 85.334407;
    c(3).constraints = 0.0071317*ants(2,:).*ants(5,:) + 0.0021813*ants(3,:).*ants(3,:) + 0.0029955*ants(1,:).*ants(2,:) - 29.48751;
    c(4).constraints = -0.0071317*ants(2,:).*ants(5,:) - 0.0021813*ants(3,:).*ants(3,:) - 0.0029955*ants(1,:).*ants(2,:) + 9.48751;
    c(5).constraints = 0.0047026*ants(3,:).*ants(5,:) + 0.0019085*ants(3,:).*ants(4,:) + 0.0012547*ants(1,:).*ants(3,:) - 15.699039;
    c(6).constraints = -0.004702*ants(3,:).*ants(5,:) - 0.0019085*ants(3,:).*ants(4,:) - 0.0012547*ants(1,:).*ants(3,:) + 10.699039;

    f = 37.293239*ants(1,:) + 0.8356891*ants(1,:).*ants(5,:)  + 5.3578547*ants(3,:).*ants(3,:) - 40792.141;
    
    %Code to calculate the residual function based on l2 norm and the
    %oracle parameter alpa which is decided on the basis of ech solution
    %constructed by the ants
    res = zeros(1,number_of_ants);
    for i = 1:number_of_ants
        res = 0;
        for j = 1:size(c,2)
            penalty(j) = min(0,c(j).constraints(i))^2;                      %l2 norm used, l1 norm may also be used here
        end
        res(1,i) = sqrt(sum(penalty));                                      %code to calculate the residual function
        if (f(i)>omega || res(1,i)>0)
            if (f(i)>omega && res(1,i)<(abs(f(i)-omega)/3))
                alpha1 = ((abs(f(i)-omega)*(1-(1/(3*sqrt(3)))))-res(1,i))/(abs(f(i)-omega)-(res(1,i)));
            elseif (f(i)>omega && res(1,i)>(abs(f(i)-omega)/3) && res(1,i)<abs(f(i)-omega))
                alpha1 = 1 - (1/(2*sqrt(abs(f(i)-omega)/res(1,i))));
            elseif (f(i)>omega && res(1,i)>abs(f(i)-omega))
                alpha1 = (1/(2*sqrt(abs(f(i)-omega)/res(1,i))));
            end
            p(1,i) = (alpha1*abs(f(i)-omega)) + ((1-alpha1)*res(1,i));      %Definition of the penalty function
        elseif (f(i)<=omega && res(1,i)==0)
            p(1,i) = -abs(f(i)-omega);
        end
    end


    best1 = [];best2 = [];best3 = [];best4 = [];best5 = [];
    [pbest1,I] = min(p);
    best1 = [best1 find( x(1,:) == ants(1,I))];                               
    best2 = [best2 find( x(2,:) == ants(2,I))];
    best3 = [best3 find( x(3,:) == ants(3,I))];
    best4 = [best4 find( x(4,:) == ants(4,I))];
    best5 = [best5 find( x(5,:) == ants(5,I))];
    pbest = [best1;best2;best3;best4;best5];
    pworst = max(p);
    
   
    for i = 1:number_of_variables
    newtau(i,:) = (1-rho).*tau(i,:);
        %for j =1:length(pbest(i,:))
            newtau(i,pbest(i,1)) = tau(i,pbest(i,1)) + eta*pbest1/pworst;  
        %end
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
        
    %Code to plot each iteration
    for i=1:number_of_ants
        plot(ants(1,i),ants(2,i),'og');
        hold on;
    end
    hold off;
    pause
  
    end
    omega = pbest;
    end
    