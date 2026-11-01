function SemiSpec_Toolkit()
%SEMISPEC_TOOLKIT  Spectral and Laplacian topological indices of semigraphs.
%
%   Companion code for the article
%   "SemiSpec: Eigenvalue- and Laplacian-Based Topological Indices of
%    Semigraphs, Resistance-Distance Formulas and Robustness of
%    Semigraph Network Models"
%   S. Naik, Author 2, Author 3, Author 4.
%
%   Running SemiSpec_Toolkit (MATLAB R2016b+ or GNU Octave 6+) will
%     1. build the adjacency graph G_a(S) and consecutive-adjacency graph
%        G_c(S) of any semigraph S given as a list of ordered s-edges;
%     2. compute 17 indices: n, m, spectral radius, energy, Estrada index,
%        natural connectivity, Laplacian energy, signless Laplacian energy,
%        Laplacian Estrada index, algebraic connectivity, Laplacian spectral
%        radius, Kirchhoff index, number of spanning trees, M1, M2, Randic
%        and sum-connectivity indices;
%     3. verify every closed-form theorem of the paper against the
%        numerical spectra (prints maximum errors);
%     4. reproduce Tables 2-6 (CSV files) and Figures 2-7 (PNG files)
%        in the folder ./SemiSpec_out .
%
%   A semigraph is entered as   p      : number of vertices (1..p)
%                               edges  : cell array of row vectors, each an
%                                        ordered s-edge (e1,e2,...,er).
%   Example (Example 3.1 of the paper):
%       p = 6; edges = {[1 2 3 4],[4 5 6],[2 6]};
%       [Aa,Ac] = sg_build(p,edges); T = sg_indices(Aa)
%
%   No toolboxes are required.

clc; close all;
outdir = fullfile(pwd,'SemiSpec_out');
if ~exist(outdir,'dir'), mkdir(outdir); end
fprintf('==============================================================\n');
fprintf(' SemiSpec toolkit : spectral indices of semigraphs\n');
fprintf('==============================================================\n\n');

%% ---------------- 1. Worked example ------------------------------------
p = 6; edges = {[1 2 3 4],[4 5 6],[2 6]};
check_semigraph(p,edges);
[Aa,Ac] = sg_build(p,edges);
Ta = sg_indices(Aa); Tc = sg_indices(Ac);
fprintf('Example 3.1  S = {(1,2,3,4),(4,5,6),(2,6)}\n');
fprintf('  adjacent degrees d_a       : %s\n', mat2str(sum(Aa,2)'));
fprintf('  consecutive degrees d_c    : %s\n', mat2str(sum(Ac,2)'));
fprintf('  A(G_a) spectrum            : %s\n', mat2str(round(1e4*Ta.lam')/1e4));
fprintf('  L(G_a) spectrum            : %s\n', mat2str(round(1e4*Ta.mu')/1e4));
print_indices('G_a(S)',Ta); print_indices('G_c(S)',Tc);
write_example_csv(fullfile(outdir,'Table2_example.csv'),Ta,Tc);

%% ---------------- 2. Complete semigraphs E_r^c and T^1_{r-1} -----------
fprintf('\nComplete semigraphs: closed form vs. numerics\n');
rs = 3:20; err = zeros(1,8); C = zeros(numel(rs),19);
for t = 1:numel(rs)
    r = rs(t);
    [pE,eE] = fam_Erc(r);  [AaE,AcE] = sg_build(pE,eE);
    [pT,eT] = fam_T1(r);   [AaT,AcT] = sg_build(pT,eT);
    IE = sg_indices(AcE); IT = sg_indices(AcT); IK = sg_indices(AaE);
    % Theorem 4.1 (G_a = K_r)
    err(1) = max(err(1), abs(IK.E  - 2*(r-1)));
    err(2) = max(err(2), abs(IK.Kf - (r-1)));
    err(3) = max(err(3), abs(IK.LEE - (1+(r-1)*exp(r)))/IK.LEE);
    err(4) = max(err(4), norm(AaE-AaT,'fro'));           % G_a(E)=G_a(T)=K_r
    % Theorem 4.2 (G_c = P_r)
    err(5) = max(err(5), abs(IE.E  - energy_path(r)));
    err(6) = max(err(6), abs(IE.EE - estrada_path(r)));
    % Theorem 4.3 (G_c = fan F_{r-1})
    err(7) = max(err(7), abs(IT.tau - fibo(2*r-2))/fibo(2*r-2));
    err(8) = max(err(8), abs(IT.Kf  - kirchhoff_fan(r)));
    C(t,:) = [r IE.E IT.E IK.E IE.EE IT.EE IK.EE IE.LE IT.LE IE.Kf IT.Kf IK.Kf ...
              IE.tau IT.tau IK.tau IE.a IT.a IK.a IT.rho];
end
fprintf('  max errors: E(K_r)=%.1e Kf(K_r)=%.1e LEE(K_r)=%.1e  Ga(E)=Ga(T): %.1e\n',err(1:4));
fprintf('              E(P_r)=%.1e EE(P_r)=%.1e tau(fan)=%.1e Kf(fan)=%.1e\n',err(5:8));
hdr = {'r','E_Gc_E','E_Gc_T','E_Kr','EE_Gc_E','EE_Gc_T','EE_Kr','LE_Gc_E','LE_Gc_T', ...
       'Kf_Gc_E','Kf_Gc_T','Kf_Kr','tau_Gc_E','tau_Gc_T','tau_Kr','a_Gc_E','a_Gc_T','a_Kr','rho_Gc_T'};
write_csv(fullfile(outdir,'Table3_complete.csv'),hdr,C);

figure('Position',[50 50 1100 580],'Color','w');
lab = {'(a) Energy','(b) Estrada index','(c) Laplacian energy', ...
       '(d) Kirchhoff index','(e) Spanning trees','(f) Algebraic connectivity'};
cols = {[2 3 4],[5 6 7],[8 9 0],[10 11 12],[13 14 15],[16 17 18]};
lg   = [0 1 0 1 1 1];
for s = 1:6
    subplot(2,3,s); hold on; box on; grid on;
    c = cols{s};
    plot(rs,C(:,c(1)),'o-','Color',[.12 .35 .65],'MarkerSize',4,'LineWidth',1.3);
    plot(rs,C(:,c(2)),'s-','Color',[.78 .33 .24],'MarkerSize',4,'LineWidth',1.3);
    if c(3) > 0, plot(rs,C(:,c(3)),'^--','Color',[.33 .33 .33],'MarkerSize',4,'LineWidth',1.3); end
    if lg(s), set(gca,'YScale','log'); end
    title(lab{s}); xlabel('r');
    if s == 1, legend({'G_c(E^c_r)=P_r','G_c(T^1_{r-1})=F_{r-1}','G_a=K_r'},'Location','northwest'); end
end
save_fig(fullfile(outdir,'Fig2_complete.png'));

%% ---------------- 3. Network families: s-chain, s-hub, s-ring ----------
fprintf('\nNetwork families (k = 4): closed form vs. numerics\n');
k = 4; ns = 2:20; fams = {'chain','hub','ring'};
NET = cell(1,3); e2 = zeros(1,7);
for f = 1:3
    M = nan(numel(ns),11);
    for t = 1:numel(ns)
        n = ns(t);
        if f == 3 && n < 3, continue; end
        [pp,ee] = family(fams{f},n,k);
        [A1,A2] = sg_build(pp,ee); I = sg_indices(A1); J = sg_indices(A2);
        M(t,:) = [n pp I.Kf I.Kf/nchoosek(pp,2) I.a I.NC I.E I.LE I.tau J.Kf I.rho];
        switch fams{f}
            case 'chain'
                e2(1) = max(e2(1), abs(I.Kf - kf_chain(n,k)));
                e2(2) = max(e2(2), abs(I.tau - (k^(k-2))^n)/I.tau);
                e2(7) = max(e2(7), abs(I.Kf - kf_block(pp,ee)));
            case 'hub'
                H = hub_closed(n,k);
                e2(3) = max(e2(3), max(abs([I.E-H.E, I.Kf-H.Kf, I.LE-H.LE, I.rho-H.rho])));
                e2(4) = max(e2(4), abs(J.Kf - kf_spider(n,k-1)));
                e2(7) = max(e2(7), abs(I.Kf - kf_block(pp,ee)));
            case 'ring'
                e2(5) = max(e2(5), abs(I.Kf - kf_ring(n,k)));
                e2(6) = max(e2(6), abs(I.tau - n*2*k^(k-3)*k^((k-2)*(n-1)))/I.tau);
        end
    end
    NET{f} = M;
    write_csv(fullfile(outdir,sprintf('Table5_%s_k4.csv',fams{f})), ...
        {'n','p','Kf','Rbar','a','NC','E','LE','tau','Kf_Gc','rho'},M);
end
fprintf('  Kf(chain)=%.1e tau(chain)=%.1e hub(E,Kf,LE,rho)=%.1e Kf(spider)=%.1e\n',e2(1:4));
fprintf('  Kf(ring)=%.1e tau(ring)=%.1e  block-resistance formula=%.1e\n',e2(5:7));

figure('Position',[50 50 1100 580],'Color','w');
ttl = {'(a) Kf(G_a)','(b) mean resistance','(c) algebraic connectivity', ...
       '(d) natural connectivity','(e) energy','(f) Kf(G_c)'};
idx = [3 4 5 6 7 10]; lgs = [1 0 1 0 0 1];
clr = {[.12 .35 .65],[.78 .33 .24],[.16 .62 .44]}; mk = {'o-','s-','^-'};
for s = 1:6
    subplot(2,3,s); hold on; box on; grid on;
    for f = 1:3
        M = NET{f}; plot(M(:,1),M(:,idx(s)),mk{f},'Color',clr{f},'MarkerSize',4,'LineWidth',1.3);
    end
    if lgs(s), set(gca,'YScale','log'); end
    title(ttl{s}); xlabel('number of s-edges n');
    if s == 1, legend({'s-chain','s-hub','s-ring'},'Location','northwest'); end
end
save_fig(fullfile(outdir,'Fig3_networks.png'));

%% ---------------- 4. Effect of cardinality k (n = 6) -------------------
ks = 2:10; figure('Position',[50 50 1100 300],'Color','w');
K = cell(1,3);
for f = 1:3
    M = zeros(numel(ks),5);
    for t = 1:numel(ks)
        [pp,ee] = family(fams{f},6,ks(t)); I = sg_indices(sg_build(pp,ee));
        M(t,:) = [ks(t) pp I.Kf/nchoosek(pp,2) I.NC I.LE/pp];
    end
    K{f} = M;
end
tt = {'(a) mean resistance','(b) natural connectivity','(c) LE(G_a)/p'};
for s = 1:3
    subplot(1,3,s); hold on; box on; grid on;
    for f = 1:3, plot(ks,K{f}(:,s+2),mk{f},'Color',clr{f},'MarkerSize',4,'LineWidth',1.3); end
    title(tt{s}); xlabel('s-edge cardinality k (n = 6)');
    if s == 1, legend({'s-chain','s-hub','s-ring'}); end
end
save_fig(fullfile(outdir,'Fig4_cardinality.png'));

%% ---------------- 5. Spectra --------------------------------------------
figure('Position',[50 50 1100 300],'Color','w');
for s = 1:2
    subplot(1,2,s); hold on; box on; grid on;
    for f = 1:3
        [pp,ee] = family(fams{f},6,4); I = sg_indices(sg_build(pp,ee));
        if s == 1, v = I.lam; else, v = I.mu; end
        plot(v, f + 0.24*(rand(size(v))-0.5),'o','Color',clr{f},'MarkerFaceColor',clr{f},'MarkerSize',4);
    end
    set(gca,'YTick',1:3,'YTickLabel',{'C_6(4)','H_6(4)','R_6(4)'}); ylim([0.5 3.5]);
    if s == 1, title('(a) adjacency spectrum of G_a'); else, title('(b) Laplacian spectrum of G_a'); end
    xlabel('eigenvalue');
end
save_fig(fullfile(outdir,'Fig5_spectra.png'));

%% ---------------- 6. Random ensemble ------------------------------------
fprintf('\nRandom ensemble (200 s-trees + 200 s-cyclic semigraphs)\n');
try, rng(2026); catch, rand('seed',2026); end %#ok<RAND>
NS = 400; X = zeros(NS,11); B = zeros(NS,9); maxerr = 0;
for t = 1:NS
    [pp,ee] = random_semigraph(t > NS/2);
    [A1,~] = sg_build(pp,ee); I = sg_indices(A1);
    X(t,:) = [I.M1 I.M2 I.R I.SCI I.E I.LE log(I.EE) log(I.LEE) I.Kf I.NC I.a];
    d = sum(A1,2); [ii,jj] = find(triu(A1)); kmax = max(cellfun(@numel,ee));
    B(t,:) = [I.n I.m I.E abs(det(A1)) I.rho kmax I.mu1 max(d(ii)+d(jj)) max(d)];
    if t <= NS/2, maxerr = max(maxerr, abs(I.Kf - kf_block(pp,ee))); end
end
R = corrcoef(X);
lbl = {'M1','M2','R','SCI','E','LE','lnEE','lnLEE','Kf','NC','a'};
write_csv(fullfile(outdir,'Table6_correlations.csv'),[{'index'} lbl],[(1:11)' R]);
fprintf('  block-resistance Kf formula on s-trees: max error %.2e\n',maxerr);
viol = sum(B(:,3) > sqrt(2*B(:,1).*B(:,2))+1e-8) + sum(B(:,3) < 2*sqrt(B(:,2))-1e-8) + ...
       sum(B(:,3) < sqrt(2*B(:,2)+B(:,1).*(B(:,1)-1).*B(:,4).^(2./B(:,1)))-1e-8) + ...
       sum(B(:,5) < B(:,6)-1-1e-8) + sum(B(:,7) < B(:,9)+1-1e-8) + sum(B(:,7) > B(:,8)+1e-8);
fprintf('  violations of Theorem 5.3 / 5.4 bounds: %d\n',viol);
c = polyfit(X(:,3),X(:,5),1); yh = polyval(c,X(:,3));
r2 = 1 - sum((X(:,5)-yh).^2)/sum((X(:,5)-mean(X(:,5))).^2);
fprintf('  energy ~ Randic : E = %.4f R + %.4f,  R^2 = %.4f\n',c(1),c(2),r2);

figure('Position',[50 50 1100 420],'Color','w');
subplot(1,2,1); imagesc(R,[-1 1]); colormap(redblue(64)); colorbar; axis square;
set(gca,'XTick',1:11,'XTickLabel',lbl,'YTick',1:11,'YTickLabel',lbl);
for i = 1:11, for j = 1:11
    text(j,i,sprintf('%.2f',R(i,j)),'HorizontalAlignment','center','FontSize',6);
end, end
title('(a) Pearson correlations');
subplot(1,2,2); hold on; box on; grid on;
plot(X(1:NS/2,3),X(1:NS/2,5),'o','Color',clr{1},'MarkerSize',3);
plot(X(NS/2+1:end,3),X(NS/2+1:end,5),'^','Color',clr{2},'MarkerSize',3);
xx = linspace(min(X(:,3)),max(X(:,3)),50); plot(xx,polyval(c,xx),'k--');
xlabel('Randic index R(G_a)'); ylabel('energy E(G_a)');
legend({'s-trees','s-cyclic',sprintf('fit, R^2 = %.3f',r2)},'Location','northwest');
title('(b) energy vs Randic index');
save_fig(fullfile(outdir,'Fig6_ensemble.png'));

figure('Position',[50 50 1100 320],'Color','w');
[~,o] = sort(B(:,3)); subplot(1,2,1); hold on; box on; grid on;
plot(sqrt(2*B(o,1).*B(o,2)),'Color',clr{2}); plot(B(o,3),'k','LineWidth',1.5); plot(2*sqrt(B(o,2)),'Color',clr{1});
legend({'sqrt(2 p m_a)','E(G_a)','2 sqrt(m_a)'},'Location','northwest'); title('(a) energy bounds');
xlabel('semigraph (sorted by energy)');
[~,o] = sort(B(:,7)); subplot(1,2,2); hold on; box on; grid on;
plot(B(o,8),'Color',clr{2}); plot(B(o,7),'k','LineWidth',1.5); plot(B(o,9)+1,'Color',clr{1});
legend({'max(d_u+d_v)','\mu_1(G_a)','\Delta_a+1'},'Location','northwest'); title('(b) Laplacian spectral radius bounds');
xlabel('semigraph (sorted by \mu_1)');
save_fig(fullfile(outdir,'Fig7_bounds.png'));

fprintf('\nAll outputs written to %s\n',outdir);
end

%% =======================================================================
%%                         CORE FUNCTIONS
%% =======================================================================
function [Aa,Ac] = sg_build(p,edges)
%SG_BUILD adjacency graph G_a (s-edges -> cliques) and consecutive graph G_c.
Aa = zeros(p); Ac = zeros(p);
for t = 1:numel(edges)
    e = edges{t};
    Aa(e,e) = 1;
    for i = 1:numel(e)-1
        Ac(e(i),e(i+1)) = 1; Ac(e(i+1),e(i)) = 1;
    end
end
Aa(1:p+1:end) = 0;
end

function check_semigraph(p,edges)
%CHECK_SEMIGRAPH verifies Sampathkumar's axioms: distinct vertices in an
%edge, |e|>=2, and two s-edges share at most one vertex.
for t = 1:numel(edges)
    e = edges{t};
    assert(numel(e) >= 2 && numel(unique(e)) == numel(e) && all(e >= 1 & e <= p), ...
        'edge %d is not a valid s-edge',t);
    for u = t+1:numel(edges)
        assert(numel(intersect(e,edges{u})) <= 1,'edges %d and %d share >1 vertex',t,u);
    end
end
end

function T = sg_indices(A)
%SG_INDICES spectral and degree-based indices of the graph with adjacency A.
n = size(A,1); d = sum(A,2); m = sum(d)/2;
lam = sort(eig((A+A')/2),'descend');
L = diag(d) - A; Q = diag(d) + A;
mu = sort(eig((L+L')/2),'descend'); mu(abs(mu) < 1e-9) = 0;
q  = eig((Q+Q')/2);
dbar = 2*m/n; nz = mu(mu > 1e-9);
T.n = n; T.m = m; T.lam = lam; T.mu = mu;
T.rho = lam(1);
T.E   = sum(abs(lam));
T.EE  = sum(exp(lam));
T.NC  = log(T.EE/n);
T.LE  = sum(abs(mu - dbar));
T.SLE = sum(abs(q - dbar));
T.LEE = sum(exp(mu));
T.a   = mu(end-1);
T.mu1 = mu(1);
if numel(nz) == n-1
    T.Kf = n*sum(1./nz); T.tau = exp(sum(log(nz)) - log(n));
else
    T.Kf = Inf; T.tau = 0;
end
[i,j] = find(triu(A));
T.M1  = sum(d.^2);
T.M2  = sum(d(i).*d(j));
T.R   = sum((d(i).*d(j)).^(-1/2));
T.SCI = sum((d(i)+d(j)).^(-1/2));
end

function print_indices(name,T)
fprintf('  %-7s n=%d m=%d rho=%.4f E=%.4f EE=%.4f NC=%.4f LE=%.4f SLE=%.4f\n', ...
    name,T.n,T.m,T.rho,T.E,T.EE,T.NC,T.LE,T.SLE);
fprintf('          LEE=%.4f a=%.4f mu1=%.4f Kf=%.4f tau=%g M1=%d M2=%d R=%.4f SCI=%.4f\n', ...
    T.LEE,T.a,T.mu1,T.Kf,round(T.tau),T.M1,T.M2,T.R,T.SCI);
end

%% ---------------- semigraph families ----------------------------------
function [p,e] = fam_Erc(r)
p = r; e = {1:r};
end

function [p,e] = fam_T1(r)
p = r; e = cell(1,r); e{1} = 1:r-1;
for i = 1:r-1, e{i+1} = [i r]; end
end

function [p,e] = family(name,n,k)
switch name
    case 'chain'          % n s-edges of cardinality k in a line
        e = cell(1,n); v = 1;
        for i = 1:n, e{i} = v:v+k-1; v = v+k-1; end
        p = v;
    case 'hub'            % n s-edges of cardinality k sharing end vertex 1
        e = cell(1,n); v = 2;
        for i = 1:n, e{i} = [1 v:v+k-2]; v = v+k-1; end
        p = v-1;
    case 'ring'           % n s-edges of cardinality k closed in a cycle
        p = n*(k-1); e = cell(1,n);
        for i = 1:n, s = (i-1)*(k-1); e{i} = mod(s+(0:k-1),p)+1; end
end
end

function [p,e] = random_semigraph(cyclic)
k = randi([2 5]); e = {1:k}; p = k; target = randi([10 40]);
while p < target
    k = randi([2 5]); anchor = randi(p);
    ed = [anchor, p+1:p+k-1]; p = p+k-1;
    if rand < 0.5, ed = fliplr(ed); end
    e{end+1} = ed; %#ok<AGROW>
end
if cyclic
    for c = 1:randi([1 3])
        Aa = sg_build(p,e);
        for tries = 1:50
            uv = randperm(p,2);
            if Aa(uv(1),uv(2)) == 0, break; end
        end
        if Aa(uv(1),uv(2)) ~= 0, continue; end
        mm = randi([0 2]); e{end+1} = [uv(1), p+1:p+mm, uv(2)]; p = p+mm; %#ok<AGROW>
    end
end
end

%% ---------------- closed forms (Section 4) ----------------------------
function E = energy_path(n)
if mod(n,2) == 0, E = 2/sin(pi/(2*(n+1))) - 2; else, E = 2*cot(pi/(2*(n+1))) - 2; end
end

function EE = estrada_path(n)
s = besseli(0,2);
for j = 1:5, s = s + 2*besseli(2*j*(n+1),2); end
EE = (n+1)*s - cosh(2);
end

function f = fibo(n)
a = 0; b = 1; for i = 1:n, [a,b] = deal(b,a+b); end; f = a;
end

function l = lucas(n)
a = 2; b = 1; for i = 1:n, [a,b] = deal(b,a+b); end; l = a;
end

function K = kirchhoff_fan(r)
n = r-1; K = 1 + r*n*lucas(2*n)/(5*fibo(2*n)) - 3*r/5;
end

function K = kf_chain(n,k)
K = n*(k-1)*((k-1)*n^2 + 3*(k-1)*n - k + 4)/(3*k);
end

function K = kf_ring(n,k)
K = n*((k-1)^2*n^2 + 6*(k-1)*(k-2)*n - (k^2+k-5))/(6*k);
end

function K = kf_spider(n,l)
K = n*(l^3-l)/6 + n*l*(l+1)/2 + n*(n-1)/2*l^2*(l+1);
end

function H = hub_closed(n,k)
N = 1+n*(k-1); D = (k-2)^2 + 4*n*(k-1);
H.E   = (2*n-1)*(k-2) + sqrt(D);
H.rho = ((k-2)+sqrt(D))/2;
H.Kf  = N*(n-1) + N*n*(k-2)/k + 1;
db    = n*k*(k-1)/N;
H.LE  = db + (n-1)*(db-1) + n*(k-2)*(k-db) + (N-db);
end

function K = kf_block(p,edges)
%KF_BLOCK Theorem 5.5: Kf(G_a) of an s-tree from branch sizes.
Aa = sg_build(p,edges); K = 0;
for t = 1:numel(edges)
    e = edges{t}; G = Aa; G(e,e) = 0;
    sz = zeros(1,numel(e));
    for i = 1:numel(e), sz(i) = numel(component(G,e(i))); end
    K = K + (2/numel(e))*(sum(sz)^2 - sum(sz.^2))/2;
end
end

function c = component(A,s)
seen = false(1,size(A,1)); seen(s) = true; stack = s;
while ~isempty(stack)
    v = stack(end); stack(end) = [];
    nb = find(A(v,:) & ~seen); seen(nb) = true; stack = [stack nb]; %#ok<AGROW>
end
c = find(seen);
end

%% ---------------- utilities -------------------------------------------
function write_csv(fname,hdr,M)
fid = fopen(fname,'w'); fprintf(fid,'%s,',hdr{1:end-1}); fprintf(fid,'%s\n',hdr{end});
for i = 1:size(M,1), fprintf(fid,'%.10g,',M(i,1:end-1)); fprintf(fid,'%.10g\n',M(i,end)); end
fclose(fid);
end

function write_example_csv(fname,Ta,Tc)
f = {'n','m','rho','E','EE','NC','LE','SLE','LEE','a','mu1','Kf','tau','M1','M2','R','SCI'};
fid = fopen(fname,'w'); fprintf(fid,'index,G_a,G_c\n');
for i = 1:numel(f), fprintf(fid,'%s,%.6f,%.6f\n',f{i},Ta.(f{i}),Tc.(f{i})); end
fclose(fid);
end

function save_fig(fname)
drawnow;
if exist('exportgraphics','file') == 2 || exist('exportgraphics','builtin') == 5
    exportgraphics(gcf,fname,'Resolution',200);
else
    print(gcf,fname,'-dpng','-r200');
end
end

function c = redblue(m)
x = linspace(-1,1,m)'; c = zeros(m,3);
c(:,1) = min(1, 1 + x).*(x < 0) + 1*(x >= 0);
c(:,3) = 1*(x <= 0) + min(1, 1 - x).*(x > 0);
c(:,2) = 1 - abs(x);
c(x < 0,1) = 1 + x(x < 0); c(x > 0,3) = 1 - x(x > 0);
c(x < 0,3) = 1; c(x > 0,1) = 1;
end
