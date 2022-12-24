%filename = {'filter_resonances_no_NOA_r9970nm_w700nm.mat','filter_resonances_NOA_index_1p56_r9970nm_w700nm.mat','filter_resonances_NOA_index_1p65_r9970nm_w700nm.mat','filter_resonances_NOA_index_1p70_r9970nm_w700nm.mat'};
filename = {'filter_resonances_no_NOA_r9970nm_w700nm.mat','filter_resonances_NOA_index_1p315_r9970nm_w700nm.mat','filter_resonances_NOA_index_1p48_r9970nm_w700nm.mat','filter_resonances_NOA_index_1p70_r9970nm_w700nm.mat'};
for ii=1:length(filename) 
    load(filename{ii});
    y = ones(1,length(lambda_res));
    stem(lambda_res,y,'LineWidth',2);
    %legend(filename{ii});
%     for kk = 1:length(lambda_res)
%         text(lambda_res(kk),y,num2str(lambda_res(kk)));
%     end
    hold on
end
legend('no NOA','NOA index 1.315','NOA index 1.48','NOA index 1.7')
title('ROut=9970nm, Width=700nm')
set(gcf,'Color','w')
%imagesc(SZ.N.x, SZ.N.y, (SZ.N.n).'); set(gca, 'YDir', 'normal');
%colorbar; %creates a picture of the layer stack
