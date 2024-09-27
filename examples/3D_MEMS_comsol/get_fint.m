function [fnl,dfnl] = get_fint(K,model,Null, Nullf, ud)

mphmesh(model)

fnl = @(input) -internalforce(model,input,Null, Nullf, ud) - K*input; % 
% fnl = @(input) 0*input;
dfnl = @(input) internalforceJacobian(model,input,Null,Nullf,ud) - K ;
end

