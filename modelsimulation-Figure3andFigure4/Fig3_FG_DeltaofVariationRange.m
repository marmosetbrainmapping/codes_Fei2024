
%%
lines = ["-","--",":","-.","-","--"];
markers = ["o","^","diamond","*",">","<"];
colors = ["r","g","b","c","y","m"];
sz = 80;

vr1 = load('HeteroVariationRange.mat').variation_range;
vr1 = squeeze(mean(vr1,1));
vr2 = load('HomoVariationRange.mat').variation_range;
vr2 = squeeze(mean(vr2,1));
vr = vr1-vr2;

figure
hold on
legends = [];
plots = [];
for i = 1:length(noises)-1
    x = Nregion_range';
    y = vr(i,:)';
    
    f = fit(x,y,'poly2');
    a = f.p1;
    b = f.p2;
    c = f.p3;
    xi = linspace(min(x),max(x),500);
    yi = a*xi.^2 + b*xi + c;
    
    scatter(x,y,sz,markers(i),'MarkerEdgeColor',colors(i),'MarkerFaceColor','none');
    p = plot(xi,yi,colors(i),'LineWidth',1.5,'LineStyle',lines(i));
    plots = [plots,p];
    legends = [legends,string(num2str(noises(i)))];
end
l = legend(plots,legends);
set(l,'Position',[0.7448    0.1958    0.1314    0.1775]);
set(gca,'Color','none','Fontsize',15);
set(l,'Color','none');
set(gcf,'Position',[0,0,700,400]);

%%

vr1 = load('HeteroVariationRange.mat').variation_range;
vr1 = squeeze(mean(vr1,1));
vr2 = load('PcontrasVariationRange.mat').variation_range;
vr2 = (squeeze(mean(vr2,1)));
vr = vr1-vr2;

figure
hold on
legends = [];
plots = [];
for i = 1:length(noises)-1
    x = Nregion_range';
    y = vr(i,:)';
    
    f = fit(x,y,'poly2');
    a = f.p1;
    b = f.p2;
    c = f.p3;
    xi = linspace(min(x),max(x),500);
    yi = a*xi.^2 + b*xi + c;
    
    scatter(x,y,sz,markers(i),'MarkerEdgeColor',colors(i),'MarkerFaceColor','none');
    p = plot(xi,yi,colors(i),'LineWidth',1.5,'LineStyle',lines(i));
    plots = [plots,p];
    legends = [legends,string(num2str(noises(i)))];
end
l = legend(plots,legends);
set(l,'Position',[0.7448    0.1958    0.1314    0.1775]);
set(gca,'Color','none','Fontsize',15);
set(l,'Color','none');
set(gcf,'Position',[0,0,700,400]);
