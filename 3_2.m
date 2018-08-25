% 定义遗传算法参数
NIND = 40;
MAXGEN = 500;
NVAR = 20;
PRECI = 20;
GGAP = 0.9;
trace = zero(MAXGEN, 2);

%建立区域描述器
FieldD = [rep([PRECI],[1,NVAR]);rep([-512;512],[1,NVAR]);rep([1;0;1;1],[1,NVAR])];
Chrom = crtbp(NIND, NVAR*PRECI);
gen = 0;
ObjV = Objfunl(bs2v(Chrom,FieldD))
while gen<MAXGEN,
	FitnV = ranking(ObjV);
	SelCh = select('sus',Chrom,FitnV,GGAP);
	SelCh = recombin('xovsp',SelCh,0.7);
	SelCh = mut(SelCh);
	ObjVSel = objfunl(bs2v(SelCh,FieldD));
	[Chrom ObjV] = reins(Chrom,SelCh,1,1,ObjV,ObjVSel);
	gen = gen +1;
	trace(gen,1) = min(ObjV);
	trace(gen,2) = sum(ObjV)/length(ObjV);
end
plot(trace(:,1));
hold on;
plot(trace(:,2),'-.');
grid;
legend('种群均值的变化','解的变化')



















