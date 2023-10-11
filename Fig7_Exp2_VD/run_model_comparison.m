
for i=1:3
    

    a = load(['./fits/fits_ID',num2str(i),'_pset1.mat']);%fixed buffer
    b = load(['./fits/fits_ID',num2str(i),'_pset2.mat']);%no buffer

    numParams = 1;
    [~,bic1] = aicbic(-a.fval,numParams,a.N);
    
    numParams = 1;
    [~,bic2] = aicbic(-b.fval,numParams,b.N);

    BF(i) = log10(exp(-1/2*(bic1 -bic2)));
end

%%

for i=1:3
    
    a = load(['./fits/fits_ID',num2str(i),'_pset1.mat']);%fixed buffer
    b = load(['./fits/fits_ID',num2str(i),'_pset3.mat']);%two kappas
    
    numParams = 1;
    [~,bic1] = aicbic(-a.fval,numParams,a.N);
    
    numParams = 2;
    [~,bic2] = aicbic(-b.fval,numParams,b.N);

    BF(i) = log10(exp(-1/2*(bic1 - bic2)));

end

