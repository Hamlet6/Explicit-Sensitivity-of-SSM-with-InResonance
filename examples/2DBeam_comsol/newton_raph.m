function root = newton_raph(func, x0, df ,K)

    epsilon = 0.00001;
    n = 100;
    for i = 1:n
        f0 = func(x0); %Calculating the value of function at x0
        f0_der=df(x0)+K; %Calculating the value of function derivative at x0
        y = x0-f0_der\f0; % The Formula
        err = norm(y-x0);
        if err<epsilon %checking the amount of error at each iteration
            break
        end
        x0 = y;
%         disp(num2str(i/n));
    end
root = y;
fprintf('No. of Iterations : %d\n',i);

end