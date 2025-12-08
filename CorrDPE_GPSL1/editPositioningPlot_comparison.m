plot(navSolutions.LLH_error(:,4));

h = get(gca,'Children');
set(gca,'Children',[h(1) h(3) h(2)]);



legend('Least Squares (Scalar Tracking)','Corr-DPE (0.1 chips)',...
    'Corr-DPE (0.05 chips)','FontSize',12);


