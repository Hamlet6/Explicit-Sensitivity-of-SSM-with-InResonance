function [fnl,dfnl] = get_fint(K,model,Null, Nullf, ud)

mphmesh(model)

fnl = @(input) -internalforce(model,input,Null, Nullf, ud) - K*input; 
dfnl = @(input) internalforceJacobian(model,input,Null,Nullf,ud) - K ;
end

