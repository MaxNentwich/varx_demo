addpath('../src/varx/')

fig_font = 9;

% Poster settings
poster_size = false;

if poster_size
    fig_font = 26;
end

%%
load ../data/neurolib_5min_rest_model_output.mat
% y = output(:,5*fs:end)'; % cut our the firt 5 seconds
y = output'; clear output

%
%y = y-mean(y);
%y = y./std(y);
model = varx(y,2);

%%
if poster_size
    figure('Units', 'inches', 'Position', [1,1,14.4,5])
else
    figure('Units', 'inches', 'Position', [1,1,8,2.5])
end

h(1)=subplot(1,4,1); 
imagesc(sqrt(Cmat)); title('Struct. Connectivity C')
ylabel('Brain areas'); xlabel('Brain areas');
colormap('hot')
axis square

h(2)=subplot(1,4,2);
R = model.A_Rvalue-diag(diag(model.A_Rvalue));
imagesc(R);  title('VARX Estimate R'); 
xlabel('Brain areas'); 
colormap('hot')
axis square

h(3)=subplot(1,4,3);
plot(sqrt(Cmat(:)),R(:),'.'); 
if poster_size
    xlabel('True'); ylabel('VARX')
else
     xlabel('True C'); ylabel('Model R')
end

title('Validation')
axis square

h(4)=subplot(1,4,4); 
lambda = 0.01;
SC = graphicalLasso(cov(y), lambda,200); 
SC = abs(SC - diag(diag(SC))); % structural connectivity
imagesc(SC); title('Sparse Inverse FC')
xlabel('Brain areas'); ylabel('Brain areas')
colormap('hot')
axis square

%%
display(['inverse cov, Spearman r=' num2str(corr(SC(:),Cmat(:),'Type','Spearman'),2)])
display(['varx, Spearman r='        num2str(corr(R(:),Cmat(:),'Type','Spearman'),2)])

fontsize(gcf, fig_font, 'points')

exportgraphics(gcf,'../results/figures/fig2_neurolib_vs_varx.png','Resolution',600)


