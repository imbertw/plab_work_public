[lambda_res omegas gammas SZ] = ring_resonances(2.0, 18.8, 1.55, 'gf45rfsoi_xs', 2.0);
save('source_resonances_no_NOA_r18800nm_w2000nm.mat');
clear all
[lambda_res omegas gammas SZ] = ring_resonances(2.0, 18.8, 1.55, 'gf45rfsoi_xs_NOA1315', 2.0);
save('source_resonances_NOA_index_1p315_r18800nm_w2000nm.mat');
clear all
[lambda_res omegas gammas SZ] = ring_resonances(2.0, 18.8, 1.55, 'gf45rfsoi_xs_NOA138', 2.0);
save('source_resonances_NOA_index_1p38_r18800nm_w2000nm.mat');
clear all
[lambda_res omegas gammas SZ] = ring_resonances(2.0, 18.8, 1.55, 'gf45rfsoi_xs_NOA148', 2.0);
save('source_resonances_NOA_index_1p48_r18800nm_w2000nm.mat');
clear all
[lambda_res omegas gammas SZ] = ring_resonances(2.0, 18.8, 1.55, 'gf45rfsoi_xs_NOA74', 2.0);
save('source_resonances_NOA_index_1p52_r18800nm_w2000nm.mat');
clear all
[lambda_res omegas gammas SZ] = ring_resonances(2.0, 18.8, 1.55, 'gf45rfsoi_xs_NOA61', 2.0);
save('source_resonances_NOA_index_1p56_r18800nm_w2000nm.mat');
clear all
filename = {'source_resonances_no_NOA_r18800nm_w2000nm.mat','source_resonances_NOA_index_1p315_r18800nm_w2000nm.mat','source_resonances_NOA_index_1p38_r18800nm_w2000nm.mat','source_resonances_NOA_index_1p48_r18800nm_w2000nm.mat','filter_resonances_NOA_index_1p52_r18800nm_w2000nm.mat','source_resonances_NOA_index_1p56_r18800nm_w2000nm.mat'};
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
legend('no NOA','NOA index 1.315','NOA index 1.38','NOA index 1.48','NOA_index_1.52','NOA index 1.56')

title('ROut=18800nm, Width=2000nm')
set(gcf,'Color','w')
%imagesc(SZ.N.x, SZ.N.y, (SZ.N.n).'); set(gca, 'YDir', 'normal');
%colorbar; %creates a picture of the layer stack