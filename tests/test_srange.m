%
% Radial scan for superbanana plateau D11/Dp
%

%% Initialise

clear;

dorun = false;
dorunshaing = true;

testcase = 'Mt1em5';

filename = ['test_srange_n3_',testcase,'.dat'];
%filename_dqdr = ['test_srange_n3_',testcase,'_dqdr.dat'];
filename_shaing = ['D11DpShaing_n3_',testcase,'.mat'];
%smin = 2.3438E-02;
%smax = 9.9219E-01;
smin = 1.578d-5;
smax = 0.86;

%% Run script

if dorun
    outfile = fopen(['../',filename],'w');
    fprintf(outfile, '%s\tr/R0\tq\tD11/Dp\n');
    fclose(outfile);

    n = 100;
    srange = linspace(sqrt(smin), sqrt(smax), n).^2;

    for s = srange
        fprintf('s = %f\n', s)
        system(sprintf('cd .. ; build/driftorbit_test %f >> %s', s, filename));
    end

end

%% load results

f = load(filename);
%f2 = load(filename_dqdr);
%f = f(f(:,1),:); % remove outer flux surfaces which do not work
%f2 = f2(f2(:,1)<.5,:); % remove outer flux surfaces which do not work
s = f(:,1);
epsr = f(:,2);
q = f(:,3);
D11Dp = f(:,4);


%% Run Shaing script

if dorunshaing
    m_phi = 3
    a = 46
    R = 181.092
    M_t = 1e-10
    C_fac = sign(M_t)

    for l = 1:100
        try
            %invAspectRatio = epsr(7*l + 1)
            %qFactor = q(7*l + 1);
            invAspectRatio = epsr(l)
            qFactor = q(l);
            E_norm = M_t*pi*qFactor % with electric field
            %E_norm = 0 % no electric field
            %C_fac = 1

            [D11OverDpl, D12OverDpl] = calc_DOverDpl_SuperbananaPlateau_chris( ...
                m_phi,a,R,invAspectRatio,qFactor,C_fac,E_norm )

            D11DpShaing(l) =  D11OverDpl*sqrt(2);
        end
    end

    %sShaing = s(1:7:end);
    %epsrShaing = epsr(1:7:end);
    %qShaing = q(1:7:end);

    sShaing = s;
    epsrShaing = epsr;
    qShaing = q;
else
    load(filename_shaing)
end

%% Plot results

% Safety factor profile
%figure(1)
%clf
%plot(s, q)
%hold on
%plot(s, f(:,5))
%plot(s(1:end-1), diff(q)./diff(s))
%xlabel('s')
%ylabel('safety factor q')
%legend('q', 'dqds (magfie)', 'dqds (Matlab)')
%hold off
%drawnow

% D11/Dp
figure(2)
clf
plot(epsr, D11Dp)
grid on
hold on
%semilogy(epsr, f2(:,4))
plot(epsr, D11DpShaing*1d-5)
%ylim([3d5,5d7])
xlabel('1/A')
ylabel('superbanana plateau D11/Dp')
legend('Collisionless', 'Shaing')
drawnow

