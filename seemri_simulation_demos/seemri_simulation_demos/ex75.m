%%
tr=linspace(0,1000,1000);
%%
plot(tr,[2*(1-exp(-tr/1000)); (1-exp(-tr/500))])
%%
plot(tr,[1.5*(1-exp(-tr/1000)); (1-exp(-tr/500))])
%%
plot(tr,[1.5*(1-exp(-tr/1000))-(1-exp(-tr/500))]); grid on
