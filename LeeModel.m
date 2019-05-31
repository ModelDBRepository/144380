%% Jonghwan Lee's model for neuronal volume dynamics

% input argument:
% d0 = initial diameter of neuronal cell body (assumed to be a sphere) [m]
% Ie0 = stimulation current amplitude [A]
% tDur = stimulation duration [s]
% bFig = plot results? (0=no, 1=yes)

% outpu argument:
% t = time vector [s]
% Vm = membrane potential [V]
% DC = relative change in the total intracellular concentration
% DV = relative change in the cell volume

% usage example:
% The following will solve Lee's model when the stimulation current of 50 pA is applied to a cell of 8-um diameter for 1 ms
% [t Vm DC DV] = LeeModel(8e-6, 5e-11, 1e-3);


function [t Vm DC DV] = LeeModel(d0, Ie0, tDur, bFig)


%% Default input arguments

	if nargin < 4
		bFig = 1;
	end
	if nargin < 3
		tDur = 1e-3;
	end
	if nargin < 2
		Ie0 = 7e-11;
	end
	if nargin < 1
		d0 = 10e-6;
	end


%% Parameters
	
	% computation parameters
	dt = 1e-5;					% time step of computation [s]
	tStart = -10e-3;			% time to start [s]
	tEnd = 50e-3;				% time to end [s]
	tStim = 0e-3;				% time to give stimulation [s]

	t = [tStart:dt:tEnd];
	nt = length(t);

	% physical constants
	e = 1.6e-19;				% charge of a single electron [C]
	na = 6e23;					% Avogadro's number
	cw = 1e6/18;				% Concentration of water in normal condition: [mol/m^3]
	RTzF = 8.3*300/1/9.6e4;		% RT/zF
	dw = 1e-10;					% water diffusion coefficient [m^2/sec]

	% cellular parameters
	dm = 1e-9;					% thickness of the cell membrane (1 nm) [m]
	cm = 10e-9*1e6;				% capacitance of the cell membrane (10 nF/mm2) [F/m^2]
	a = 4*pi*(d0/2)^2;			% cell area
	v = 4*pi/3*(d0/2)^3;		% cell volume
	cKi = 150;	cKo = 5.5;		% intracellular/extracellular concentration of K [mol/m^3]
	cNi = 15;	cNo = 150;		% intracellular/extracellular concentration of Na [mol/m^3]
	nK = cKi*v;					% intracellular amount of K [mol]
	nN = cNi*v;					% intracellular amount of Na [mol]

	% electrophysiological parameters in the resting state
	Vm0 = -67e-3;				% membrane potential [V]
	eK0 = RTzF*log(cKo/cKi);	% equilibrium potential of K [V]
	eN0 = RTzF*log(cNo/cNi);	% equilibrium potential of Na [V]
	eL = -54.4e-3;				% equilibrium potential of leak current [V]
	gK0 = 360;					% conductivity of K ion channel [S/m^2]
	gN0 = 1200;					% conductivity of Na ion channel [S/m^2]
	gL0 = 3;					% equivalent conductivity of leak current channel [S/m^2]
	jpK0 = -0.0452;				% K current by Na-K active pump 
	jpN0 = 0.0071;				% Na current by Na-K active pump 


%% Time-varying variables
	
	% membrane potential and current
	Ie = zeros(nt,1);			% stimulation current [A]
	Vm = zeros(nt,1);			% membrane potential [V]
	JK = zeros(nt,1);			% membrane current of K ion
	JN = zeros(nt,1);			% membrane current of Na ion
	JL = zeros(nt,1);			% leak current
	JW = zeros(nt,1);			% water flux
	JPK = zeros(nt,1);			% K current by Na-K active pump
	JPN = zeros(nt,1);			% Na current by Na-K active pump

	% relative changes
	DK = zeros(nt,1);			% relative change in the intracellular amount of K ion : nK(t) = nK(t=0) * (1 + DK(t))
	DN = zeros(nt,1);			% relative change in the intracellular amount of Na ion : nN(t) = nN(t=0) * (1 + DN(t))
	DV = zeros(nt,1);			% relative change in the cell volume : v(t) = v(t=0) * (1 + DV(t))
	
	% open probabilities of subunits
	n = zeros(nt,1);
	m = zeros(nt,1);
	h = zeros(nt,1);

	% transfer rates & intracellular concentrations: no need to save to memory
	an = 0;		bn = 0;
	am = 0;		bm = 0;
	ah = 0;		bh = 0;
	cK = 0;		cN = 0;


%% Stimulation current

	for it=round((tStim-tStart)/dt):round((tStim-tStart+tDur)/dt)
		Ie(it) = Ie0;
	end


%% Solve DE using Runge-Kutta
	
	for it=1:nt

		% obtain solutions at time t from the results at time (t-1)
		if it == 1
			% membrane potential
			Vm(it) = Vm0;
			% relative changes in the intracellular ion amount
			DK(it) = 0;
			DN(it) = 0;
			% relative volume change
			DV(it) = 0;
		else
			% membrane potential
			Vm(it) = Vm(it-1) + dt/cm* ( Ie(it-1)/a/(1+2/3*DV(it-1)) - JK(it-1) - JN(it-1) - JL(it-1) );
			% relative changes in the intracellular ion amount
			DK(it) = DK(it-1) - dt/e/na/nK *(JK(it-1)+JPK(it-1)) *4*pi*(3/4/pi *v*(1+DV(it-1)))^(2/3);
			DN(it) = DN(it-1) - dt/e/na/nN *(JN(it-1)+JPN(it-1)) *4*pi*(3/4/pi *v*(1+DV(it-1)))^(2/3);
			% relative volume change
			DV(it) = DV(it-1) - dt/v/cw * JW(it-1) * a*(1+2/3*DV(it-1));
		end

		% calculate parameters at time t, which need for the calculation of membrane currents
		if it == 1
			% transfer rates
			an = 1e4 *(Vm0+0.055)/( 1-exp(-100*(Vm0+0.055)) );
			bn = 125 *exp(-12.5*(Vm0+0.065));
			am = 1e5 *(Vm0+0.040)/( 1-exp(-100*(Vm0+0.040)) );
			bm = 4000 *exp(-1000/18*(Vm0+0.065));
			ah = 70 *exp(-50*(Vm0+0.065));
			bh = 1e3 /( 1+exp(-100*(Vm0+0.035)) );
			% n, m, h
			n(it) = an/(an+bn);
			m(it) = am/(am+bm);
			h(it) = ah/(ah+bh);
			% intracellular concentration
			cK = nK/v;
			cN = nN/v;
		else
			% transfer rates
			an = 1e4 *(Vm(it-1)+0.055)/( 1-exp(-100*(Vm(it-1)+0.055)) );
			bn = 125 *exp(-12.5*(Vm(it-1)+0.065));
			am = 1e5 *(Vm(it-1)+0.040)/( 1-exp(-100*(Vm(it-1)+0.040)) );
			bm = 4000 *exp(-1000/18*(Vm(it-1)+0.065));
			ah = 70 *exp(-50*(Vm(it-1)+0.065));
			bh = 1e3 /( 1+exp(-100*(Vm(it-1)+0.035)) );
			% n, m, h
			n(it) = n(it-1) + dt* ( an*(1-n(it-1)) - bn*n(it-1) );
			m(it) = m(it-1) + dt* ( am*(1-m(it-1)) - bm*m(it-1) );
			h(it) = h(it-1) + dt* ( ah*(1-h(it-1)) - bh*h(it-1) );
			% intracellular concentration
			cK = nK*(1+DK(it-1)) / v/(1+DV(it-1));
			cN = nN*(1+DN(it-1)) / v/(1+DV(it-1));
		end

		% calculate membrane currents at time t, which will be used at the next time step (membrane currents can be directly determined, not by differential equation)
		eK = RTzF*log(cKo/cK);
		eN = RTzF*log(cNo/cN);
		JK(it) = gK0*n(it)^4 * (Vm(it)-eK);
		JN(it) = gN0*m(it)^3*h(it) * (Vm(it)-eN);
		JL(it) = gL0 * (Vm(it)-eL);
		JW(it) = -1*dw/dm* ( cK+cN -cKi-cNi);
		% on the short time scale (milliseconds), we can assume the constant active pump currents.
		JPK(it) = jpK0;
		JPN(it) = jpN0;
		% on the large time scale (minutes), we need to consider homeostasis in the intracellular concentration maintained by dynamic active pump currents.
%		tau = 10;									% time constant of homeostasis
%		JPK(it) = jpK0 + e*na*nK/a/tau*DK(it);		
%		JPN(it) = jpN0 + e*na*nN/a/tau*DN(it);

		if n(it) > 1 || n(it) < 0
			disp(['ERROR: n = ' num2str(n(it),3) ' at it = ' num2str(it)]);
			break;
		end
		if m(it) > 1 || m(it) < 0
			disp(['ERROR: m = ' num2str(m(it),3) ' at it = ' num2str(it)]);
			break;
		end
		if h(it) > 1 || h(it) < 0
			disp(['ERROR: h = ' num2str(h(it),3) ' at it = ' num2str(it)]);
			break;
		end
		
	end

	DC = ( nK*(1+DK)/v./(1+DV) + nN*(1+DN)/v./(1+DV) ) / (cKi+cNi) - 1;


%% Plot results

	if bFig == 1

		figure(1);  clf;
			subplot(2,2,1);  plot(t,Ie);  axis tight;  xlabel('t');  title('Ie');
			subplot(2,2,2);  plot(t,n.^4,'r', t*1e3,m.^3.*h,'b');  axis tight;  xlabel('t');  title('n^4 (red), m^3h (blue)');
			subplot(2,2,3);  plot(t,JK*a,'r', t,JN*a,'b', t,JL*a,'g');  axis tight;  xlabel('t');  title('I_K (red), I_{Na} (blue), I_L (green)');
			subplot(2,2,4);  plot(t,Vm,'k');  axis tight;  xlabel('t');  title('Vm');

		figure(2);  clf;
			subplot(2,2,1);  plot(t,JK-JK(1),'r', t,JN-JN(1),'b');  axis tight;  xlabel('t');  title('JK-JK(0) (red), JN-JN(0) (blue)');
			subplot(2,2,2);  plot(t,nK*(1+DK)/v./(1+DV)-cKi,'r', t,nN*(1+DN)/v./(1+DV)-cNi,'b', t,nK*(1+DK)/v./(1+DV)-cKi + nN*(1+DN)/v./(1+DV)-cNi,'m');  axis tight;  xlabel('t');  title('\Delta[K]_i (red), \Delta[Na]_i (blue), \Delta[K+Na]_i (magenta)');
			subplot(2,2,3);	 plot(t,(JPK-jpK0)/jpK0*sign(jpK0),'r', t,(JPN-jpN0)/jpN0*sign(jpN0),'b');  axis tight;  xlabel('t');  title('JPK-JPK(0) (red), JNK-JNK(0) (blue)');
			subplot(2,2,4);  plot(t,JW,'b');  axis tight;  xlabel('t');  title('J_{water}');

		figure(3);  clf;
			subplot(2,2,1);  line(t*1e3,Vm*1e3,'color','k');  axis tight;  ylim([-100 60]);  xlabel('Time [ms]');  title('Membrane Potential [mV]');
			subplot(2,2,2);  line(t*1e3,DC,'color','k');  axis tight;  ylim([-1 1]*5e-5);  xlabel('Time [ms]');  title('Relative Concentration Change');
			subplot(2,2,3);  line(t*1e3,DV,'color','k');  axis tight;  ylim([-1 1]*2e-5);  xlabel('Time [ms]');  title('Relative Volume Change');

	end

