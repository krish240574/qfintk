EU:{[s0;pu;pd;N;r;dv;iscall]
			/ European pricing
			M:1+N;  / Number of terminal nodes of tree
			u::1+pu; / Expected value in the up state
			d::1-pd; / Expected value in the down state
			qu::((exp(dt*(r-dv)))-d)%(u-d);
			qd::1-qu;
			/ Initialize stock options tree
			STs:{s0*(u xexp (N-x))*(d xexp x)}each til M
			/ Initialize payoffs tree
			payoffs:0.0;
			$[iscall=1;payoffs:0|(STs-K);payoffs:0|(K-STs)];
			/ traverse the tree backwards
			{payoffs::((payoffs[til -1+count payoffs] * qu) + payoffs[(1+til count payoffs)[til -1+count payoffs]]*qd)*df}each til N;
		};
US:{[dummy]
			/ US pricing
			u::1+pu; / Expected value in the up state
			d::1-pd; / Expected value in the down state
			qu::((exp(dt*(r-dv)))-d)%(u-d);
			qd::1-qu;
			TREE[iscall];
		};

CRR:{[dummy]
			/ CRR model for BOPM
			dt:T%N;
			u::exp(sd*sqrt(dt));
			d::1.0%u;
			qu::((exp(dt*(r-dv)))-d)%(u-d);
			qd::1-qu;
			show qu;
			show qd;
			show "CRR";
			TREE[iscall];
	};

LR:{[dummy]
	$[1=N mod 2;oddN:N;oddN:N+1];
	d1:(log(s0%K)+((r-dv)+(sd*sd)%2.0)*T)%(sd*sqrt T);
	d2:(log(s0%K)+((r-dv)-(sd*sd)%2.0)*T)%(sd*sqrt T);
	pp2inversion:{[z;n]0.5+$[z<0;-1;1]*sqrt(0.25-0.25*exp(-((z%(n+1.0%3.0+0.1%(n+1))) xexp 2 )*(n+1.0/0.6)))};
	pbar:pp2inversion[d1;oddN];
	p:pp2inversion[d2;oddN];
	u:(1%df)*(pbar*p);
	d:((1%df)-(p*u))%p;
	qu:p;
	qd:1-qu;
	TREE[iscall];
	};

TREE:{[iscall]
			/ All tree related functionality here
			show "TREE";
			STs::(1,1)#s0;

			/ Stock prices path simulation
			{pb:last STs;STs::STs,enlist (pb*u),(d*last pb)}each til N;

			/ Initialize payoffs
			payoffs::0.0;
			$[iscall=1;payoffs::0|(STs[N]-K);payoffs::0|(K-STs[N])];
			show "-------";
			show payoffs;
			show "-------";

			/ traverse tree and check for early exercise too.
			{	
				/ Payoffs from not exercising the options
				payoffs::((payoffs[til -1+count payoffs] * qu)+payoffs[(1+til count payoffs)[til -1+count payoffs]]*qd)*df;
				/ Payoffs from early exercise of options
				if[iseu=0;payoffs::payoffs|$[iscall=1;STs[x]-K;K-STs[x]]];
				show x;
				show payoffs;
			}each reverse til N;	

	};

main:{[dummy]
	s0::50f;
	pu::0.2;
	pd::0.2;
	N::2;
	r::0.05;
	dv::0.0;
	iscall::1;
	iseu::0;
	sd::0.3;
	T::0.5;
	K::50f;
	dt::T%N;
	df::exp(neg (r-dv)*dt);
	US[0];
	};

main[0];
