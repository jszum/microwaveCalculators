import math as Math


ms_c = 2.997925e8; 					# speed of light in meters / second
msmpm = 39370.07874; 				# mils per meter
ms_fsmu = 4.0 * Math.pi * 1.0e-7; 	# freespace mu
ms_rc = 5.96e7; 					#conductivity of copper.
msZofn = 59.9585;

#user vars
msEr = 10

def coth2(x):
	retVal = ((Math.exp(x) + Math.exp(-x)) / 2.0) / ((Math.exp(x) - Math.exp(-x)) / 2.0);
	return retVal * retVal;

	
def Zo1(u):

	f = 0.0;
	z = 0.0;
	f = 6.0 + ( 2.0 * Math.pi - 6.0) * Math.exp(-1.0 * Math.pow((30.666 / u), 0.7528));
	z = msZofn * Math.log(f / u + Math.pow(1.0 + 4.0 / (u * u), 0.5));
	return z;

def msCalcEeff(u, er):
	a = 0.0;
	b = 0.0;
	ef = 0.0;
	a = 1.0 + Math.log((Math.pow(u,  4.0)+Math.pow(u/52.0, 2.0)) / (Math.pow(u,  4.0) + 0.432)) / 49.0 + Math.log1p(Math.pow(u / 18.1, 3)) / 18.7;
	b = 0.564 * Math.pow((er - 0.9) / (er + 3.0), 0.053);	
	ef = (er + 1.0) / 2.0 + ((er - 1) / 2.0) * Math.pow(( 1.0 + 10.0 / u), -1.0 * a * b);
	return ef;
	

def msCalcZo(h, t, w, f):

	eff0 = 0.0

	u = w / h;
	du1 = 0;
	dur = 0;
	if (t > 0):
		tu = t / h;
		du1 = (tu / Math.pi) * Math.log1p(4.0 * Math.e / (tu * coth2(Math.pow(6.517 * u, 0.5))));
		dur = 0.5 * (1.0 + 1.0 / Math.cosh(Math.pow(msEr - 1.0, 0.5))) * du1;
	
	u1 = u + du1;
	ur = u + dur;
	Zo = Zo1(ur) / Math.pow(msCalcEeff(ur, msEr), 0.5);
	eff0 = msCalcEeff(ur, msEr) * Math.pow(Zo1(u1) / Zo1(ur), 2.0);
	fn = f * h * 0.0254 / 1.0e9;
	p1 = 0.27488 + u*(0.6315 + 0.525 / Math.pow(1.0 + 0.0157*fn,  20.0)) - 0.065683*Math.exp(-8.7513*u);
	p2 = 0.33622*(1.0 - Math.exp(-0.03442*msEr));
	p3 = 0.0363*Math.exp(-4.6*u)*(1.0 - Math.exp(-1.0*Math.pow(fn / 3.87, 4.97)));
	p4 = 2.751*(1.0 - Math.exp(-1.0*Math.pow(msEr / 15.916,  8.0))) + 1.0;
	p = p1*p2*Math.pow(fn*(0.1844 + p3*p4), 1.5763);
	msEeff = (eff0  + msEr * p) / (1.0 + p);
	r1 = 0.03891*Math.pow(msEr,  1.4);
	r2 = 0.267*Math.pow(u,  7.0);
	r3 = 4.766*Math.exp(-3.228*Math.pow(ur,  0.641));
	r4 = 0.016 + Math.pow(0.0514*msEr, 4.524);
	r5 = Math.pow(fn / 28.843,  12.0);
	r6 = 22.2 * Math.pow(ur,  1.92);
	r7 = 1.206 - 0.3144*Math.exp(-r1)*(1.0 - Math.exp(-r2));
	r8 = 1.0 + 1.275*(1.0 - Math.exp(-0.004625*r3*Math.pow(msEr, 1.674)*Math.pow(fn/18.365, 2.745)));
	r9 = (5.086*r4*r5 / (0.3838 + 0.386*r4)) * (Math.exp(-r6) / (1.0 + 1.2992*r5))*Math.pow(msEr - 1.0, 6.0) / (1.0 + 10.0*Math.pow(msEr - 1.0, 6.0));
	r10 = 0.00044 * Math.pow(msEr, 2.136) + 0.0184;
	r11 = Math.pow(fn / 19.47, 6.0) / (1.0 + 0.0962*Math.pow(fn/19.47, 6.0));
	r12 = 1.0 / (1.0 + 0.00245*u*u);
	r13 = 0.9408*Math.pow(msEeff, r8) - 0.9603;
	r14 = (0.9408 - r9)*Math.pow(eff0, r8) - 0.9603;
	r15 = 0.707*r10*Math.pow(fn/12.3, 1.097);
	r16 = 1.0 + 0.0503*msEr*msEr*r11*(1.0 - Math.exp(-1.0*Math.pow(u/15.0, 6.0)));
	r17 = r7*(1.0 - 1.1241*(r12/r16)*Math.exp(-0.026*Math.pow(fn, 1.15656)- r15));
	msZo = Zo * Math.pow(r13 / r14, r17);
	
	return msZo


print msCalcZo(1,1,1,1)