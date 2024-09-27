figure; hold on
plot_stab_lines(pFRCs{1}{1}.om,pFRCs{1}{1}.rho,pFRCs{1}{1}.st,ST1);
plot_stab_lines(pFRCs{1}{2}.om,pFRCs{1}{2}.rho,pFRCs{1}{2}.st,ST1);
plot_stab_lines(pFRCs{2}{1}.om,pFRCs{2}{1}.rho,pFRCs{2}{1}.st,ST2);
plot_stab_lines(pFRCs{2}{2}.om,pFRCs{2}{2}.rho,pFRCs{2}{2}.st,ST2);
plot_stab_lines(pFRCs{3}{1}.om,pFRCs{3}{1}.rho,pFRCs{3}{1}.st,ST3);
plot_stab_lines(pFRCs{3}{2}.om,pFRCs{3}{2}.rho,pFRCs{3}{2}.st,ST3);
xlim(freqRange); grid on; box on
xlabel('\Omega'); ylabel('\rho');

plot(frcp1{1}{1}.om(1:5:end),frcp1{1}{1}.rho(1:5:end),'kx')
plot(frcp1{1}{2}.om(1:5:end),frcp1{1}{2}.rho(1:5:end),'kx')
plot(frcp3{1}{1}.om(1:5:end),frcp3{1}{1}.rho(1:5:end),'rx')
plot(frcp3{1}{2}.om(1:5:end),frcp3{1}{2}.rho(1:5:end),'rx')


figure; hold on
plot_stab_lines(pFRCs{1}.om,pFRCs{1}.ampL2,pFRCs{1}.st,ST1);
plot_stab_lines(pFRCs{2}.om,pFRCs{2}.ampL2,pFRCs{2}.st,ST2);
plot_stab_lines(pFRCs{3}.om,pFRCs{3}.ampL2,pFRCs{3}.st,ST3);
xlim(freqRange); grid on; box on
xlabel('\Omega'); ylabel('Amp');
