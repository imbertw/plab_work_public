% [lambda_res omegas gammas SZ] = ring_resonances(0.7, 9.97, 1.55, 'gf45rfsoi_xs', 2.0);
% save('filter_resonances_no_etch_r9970nm_w700nm.mat');
% clear all
% [lambda_res omegas gammas SZ] = ring_resonances_box(0.7, 9.97, 1.55, 'gf45rfsoi_xs_box', 2.0, 0.145);
% save('filter_resonances_etch_50nm_r9970nm_w700nm.mat');
% clear all
% [lambda_res omegas gammas SZ] = ring_resonances_box(0.7, 9.97, 1.55, 'gf45rfsoi_xs_box', 2.0, 0.095);
% save('filter_resonances_etch_100nm_r9970nm_w700nm.mat');
% clear all
% [lambda_res omegas gammas SZ] = ring_resonances_box(0.7, 9.97, 1.55, 'gf45rfsoi_xs_box', 2.0, 0.045);
% save('filter_resonances_etch_150nm_r9970nm_w700nm.mat');
% clear all
%filename = {'filter_resonances_no_etch_r9970nm_w700nm.mat','filter_resonances_etch_50nm_r9970nm_w700nm.mat','filter_resonances_etch_100nm_r9970nm_w700nm.mat','filter_resonances_etch_150nm_r9970nm_w700nm.mat'};
filename = {'filter_resonances_no_etch_r9970nm_w700nm.mat','filter_resonances_etch_50nm_r9970nm_w700nm.mat','filter_resonances_etch_100nm_r9970nm_w700nm.mat', 'filter_resonances_etch_150nm_r9970nm_w700nm.mat'};
for ii=1:length(filename) 
    load(filename{ii});
    y = ones(1,length(lambda_res(gammas==70)));
    stem(lambda_res(gammas==70),y,'LineWidth',2);
    %legend(filename{ii});
%     for kk = 1:length(lambda_res)
%         text(lambda_res(kk),y,num2str(lambda_res(kk)));
%     end
    hold on
end
legend('no etch','etch 50nm', 'etch 100nm', 'etch 150nm')
title('Filter ring resonances for different oxide etch depths, ROut=9970nm, Width=700nm')
set(gcf,'Color','w')
%imagesc(SZ.N.x, SZ.N.y, (SZ.N.n).'); set(gca, 'YDir', 'normal');
%colorbar; %creates a picture of the layer stack