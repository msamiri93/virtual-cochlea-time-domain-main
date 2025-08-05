fdir = '.\houtput\';

rlist = ls(fdir);

icount = 0;
for ii = 1:size(rlist,1)
    if contains(rlist(ii,:),'r052925')
        icount = icount + 1;
        load([fdir,rlist(ii,:),'\model.mat'],'MP','Nd');
        mlist = ls([fdir,rlist(ii,:)]);
        for mm = 1:size(mlist,1)
            if contains(mlist(mm,:),'cycle')
            load([fdir,rlist(ii,:),'\',mlist(mm,:)]);
            end
        end
        tt = linspace(0,MP.period,MP.steps+1);
        uf = (2/MP.steps)*uu*exp(-1i*2*pi*MP.freq*tt');
        nAF = strcmp(Nd.name,'AF');
        nBB = strcmp(Nd.name,'BB');
        nDD = strcmp(Nd.name,'DD');        
        xx = reshape(uf,6,Nd.N);        
        R.yBM = xx(2,nAF,:);
        R.yRL = xx(2,nBB,:);
        R.xDC = xx(1,nDD,:);
        R.yDC = xx(2,nDD,:);
        if isfield(MP,'g_kA')
            R.kA = MP.g_kA;
        else
            R.kA = 1;
        end        
        rr(icount) = R;
        clear R
    end
end

%%
a_idx = find([rr.sensitive_cochlea]);
f_idx = find([rr(a_idx).freq] == 3);

[map,num] = brewermap(3+numel(f_idx),'Purples');
yy = reshape([rr(a_idx(f_idx)).yBM],1201,numel(f_idx));

% figure;
hold on
[~,s_idx] = sort([rr(a_idx(f_idx)).kA]);
for ii = 1:numel(f_idx)

    cc = map(3+s_idx(ii),:);
    if rr(a_idx(f_idx(ii))).kA == 1
        cc = [0,0,0];
    end
    hold on; plot(rr(1).xx,abs(yy(:,ii)),'Color',cc,'LineWidth',1)
end

clim([0.01,10])
colormap(map)
set(gca,'ColorScale','log')
colorbar

%%

freq = unique([rr.freq]);
a_idx = find([rr.sensitive_cochlea]);

figure;
hold on
nn = 2;
[map,num] = brewermap(nn+numel(freq),'Greens');

for ii = 1:numel(freq)
    f_idx = find([rr(a_idx).freq] == freq(ii));
    yy = reshape([rr(a_idx(f_idx)).yBM],1201,numel(f_idx));
    [yy_max] = max(abs(yy),[],1);
    kA = [rr(a_idx(f_idx)).kA];
    [~,s_idx] = sort(kA);
    plot(kA(s_idx),yy_max(s_idx),'Color',map(nn+ii,:),'LineStyle','-','Marker','o','LineWidth',1);
end

naxis
clim([1,20]);
set(gca,'ColorScale','log')
colormap(map)
colorbar
xscale('log')
yscale('log')
xlabel('x(kA)')
ylabel('y_{BM} (\mum)')