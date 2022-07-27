
x = Normal(0,1);
p = Beta(10,10);
qn = GIPT(x,p);

f1 = figure('Position',[100 100 700 300]);

xx = -3:0.01:3;
xpdf = x.PDF(xx);
sp1 = subplot(1,3,1);
plot(xx,xpdf);
xlabel('X');
ylabel('PDF(X)')
title(['(a)  X=' x.StringName],'interpreter','latex');

px = 0.2:0.01:0.8;
ppdf = p.PDF(px);
sp2 = subplot(1,3,2);
plot(px,ppdf);
xlabel('P');
ylabel('PDF(P)')
title(['(b)  P=' p.StringName],'interpreter','latex');

qx = -1:0.01:1;
qpdf = qn.PDF(qx);
sp3 = subplot(1,3,3);
plot(qx,qpdf);
% xlabel('$T=F_X^{-1}(p)$','interpreter','latex')
xlabel('T')
ylabel('PDF(T)');
title(['(c)  T=GIPT(X,P)$=F_X^{-1}(P)$'],'interpreter','latex')

set(f1,'PaperSize',[7 3]);   % Size, in inches, of the figure plotted on paper.
print(f1,'GIPTexample.pdf','-dpdf','-bestfit');