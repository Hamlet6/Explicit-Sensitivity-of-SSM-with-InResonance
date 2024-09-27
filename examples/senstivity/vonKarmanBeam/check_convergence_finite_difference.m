

% computation of peturbed SSMs
dKc = [5:-1:2 1:-0.1:0.2 0.1 0.05 0.01 0.001];
n = numel(dKc);
dw21_dm = zeros(596,n);
dw30_dm = zeros(596,n);
dlamd = zeros(n,1);
for k=1:n
    ssm_mp = comp_ssm(kLin+dKc(k),kNon);
    ssm_mm = comp_ssm(kLin-dKc(k),kNon);
    dw21_dm(:,k) = (ssm_mp.W21-ssm_mm.W21)/(2*dKc(k));
    dw30_dm(:,k) = (ssm_mp.W30-ssm_mm.W30)/(2*dKc(k));
    dlamd(k) = (ssm_mp.lamd-ssm_mm.lamd)/(2*dKc(k));
end

figure; semilogx(dKc,abs(dlamd));
figure; semilogx(dKc,abs(dw21_dm(474,:)))
figure; semilogx(dKc,abs(dw30_dm(393,:)))