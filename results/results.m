%%  PreProcessing
clear, close all
clc
%%  Initializing

ScreenSize      = get(0,'ScreenSize');
addpath('../main function/')
addpath('../functions/')
addpath('../data/')

figSelector, uiwait()

if isempty(whos('figureNum'))
    return
end

%% Figure 2 (Data Fiting)

if any(2 == figureNum)

    load SMART_Result

    minNs           = min(SMART_Result(:,2));
    maxNs           = max(SMART_Result(:,2));
    NsVec           = (minNs:maxNs).';

    DOFmin          = zeros(size(SMART_Result,2),1);

    % - Data Fiting
    counter         = 0;

    for Ns = minNs:maxNs

        counter = counter + 1;
        minIdx = find(SMART_Result(:,2) == Ns, 1 );
        DOFmin(counter) = SMART_Result(minIdx,1);

    end

    f               = fittype('a*sqrt(x) + b', 'independent', 'x', 'coefficients', {'a', 'b'});
    sqrtFit         = fit(DOFmin, NsVec, f, 'StartPoint', [minNs,DOFmin(1)]);
    fitData         = sqrtFit.a * sqrt(DOFmin) + sqrtFit.b;

    % - Plot Results
    figure('Name','Data Fiting','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/3) floor(ScreenSize(4)/4)]);
    movegui('center')
    scatter(DOFmin,NsVec,8,'filled', 'k'), hold on
    plot(DOFmin,fitData, 'LineWidth',1.5, 'color', '#0072BD')
    ylim([0 200])
    xlim([0 5e3])
    grid on
    box on
    xlabel('DOF','Interpreter','Latex','fontsize',12)
    ylabel('Number of Sensors ($N_s$)','Interpreter','Latex','fontsize',12)
    legend('Data', 'Fitted Curve','Interpreter','Latex','fontsize',12,'location','nw')

    drawnow
    clearvars -except figureNum ScreenSize

end

%%  Figure 3 (Complexity)

if any(3 == figureNum)

    load OptResult
    load SMART_Result

    % - Complexity
    compPro     = log10(OptResult(1:find(OptResult(:,1) == max(OptResult(:,1))),2).^(5));
    compMSA     = log10(OptResult(:,1).^(2*10) .* OptResult(:,2));

    % - Plot Results
    figure('Name','Complexity','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/3) floor(ScreenSize(4)/4)]);
    movegui('center')
    plot(compPro,'-*', 'LineWidth',1.5, 'MarkerIndices', 1:8:213, 'MarkerSize', 7), hold on
    plot(compMSA,'-^', 'MarkerFaceColor', '#D95319', 'LineWidth',1.5, 'MarkerIndices', 1:8:213, 'MarkerSize', 7)
    ylabel('$\log_{10}(\mathrm{complexity})$', 'Interpreter','Latex','FontSize',12)
    xlabel('DOF', 'Interpreter','Latex','FontSize',12)
    legend('(E-)SMART', 'MSA','Interpreter','Latex','fontsize',12,'location','east')
    xlim([-inf inf])
    grid minor

    drawnow
    clearvars -except figureNum ScreenSize

end

%%  Figure 4 (Convergence Curve)

if any(4 == figureNum)

    DOFvec              = [30 40 50 60];

    for i = 1:length(DOFvec)

        DOF = DOFvec(i);

        [~, objfunVal_I] = SMART(DOF, 'result', 'off');
        [~, objfunVal_II] = E_SMART(DOF, 'result', 'off');

        figure('Name',['Convergence Curve for DOF = ', num2str(DOF)],'NumberTitle','off',...
            'Position', [(i-1)*floor(ScreenSize(3)/4) floor(ScreenSize(4)-floor(ScreenSize(3)/4))/2 floor(ScreenSize(3)/4) floor(ScreenSize(3)/4)]);
        plot(0:length(objfunVal_I)-1, objfunVal_I, '-*', 'LineWidth', 1.5), hold on
        plot(0:length(objfunVal_II)-1, objfunVal_II, '-v', 'LineWidth', 1.5)
        grid on
        axis square
        box on
        legend('SMART', 'E-SMART', 'Interpreter','latex','Location','ne','FontSize',13)
        xlim([0 max(length(objfunVal_I), length(objfunVal_II))-1])
        xticks(0:max(length(objfunVal_I), length(objfunVal_II))-1)
        ylim([0 100*ceil(sum(0:DOF)/100)])
        xlabel('Inter CD Iteration','Interpreter','latex','FontSize',13)
        ylabel('Sum of Uncovered Holes','Interpreter','latex','FontSize',13)
        ax = gca;
        ax.FontSize = 12;

    end

    drawnow
    clearvars -except figureNum ScreenSize

end

%%  Figure 5 (Sparse Array Design Results)

if any(5 == figureNum)

    load SMART_Result
    load E_SMART_Result
    load OptResult
    load CCAResult

    nMax                = 32;
    nVec                = 3:nMax;

    DOFmax              = 32;

    NsSMART             = zeros(DOFmax-1, 1);
    NsE_SMART           = zeros(DOFmax-1, 1);
    DOFaug              = zeros(DOFmax-1, 1);
    NsAug               = zeros(DOFmax-1, 1);
    DOFmisc             = zeros(DOFmax-1, 1);
    NsMISC              = zeros(DOFmax-1, 1);

    % - Main Loop
    % SMART and E-SMART
    counter             = 0;
    for DOF = 2:DOFmax

        counter = counter + 1;

        indicesSMART = SMART(DOF,'result','off');
        NsSMART(counter) = length(indicesSMART);
        if DOF >= 5
            indicesE_SMART = E_SMART(DOF, 'result', 'off');
            NsE_SMART(counter) = length(indicesE_SMART);
        end

    end

    % ANA-II1
    counter             = 0;
    for nAug = min(NsSMART):max(NsSMART)

        counter = counter + 1;
        indicesAug = AugNestedII1Gen(nAug);

        if max(indicesAug) > DOFmax
            break
        end

        DOFaug(counter) = max(indicesAug);
        NsAug(counter) = length(indicesAug);

    end

    % MISC
    counter             = 0;
    for nMISC = 5:max(NsSMART)

        counter = counter + 1;
        indicesMISC = MISC(nMISC);

        if max(indicesMISC) > DOFmax
            break
        end

        DOFmisc(counter) = max(indicesMISC);
        NsMISC(counter) = length(indicesMISC);

    end

    % Complementary Co-prime Array (CCA)
    CCAResult_prime_free = CCAResult;
    CCAResult_prime_free(isprime(CCAResult_prime_free(:,1)), :) = [];

    DOFmaxAug           = zeros(size(nVec));
    DOFmaxMISC          = zeros(size(nVec));

    for n =  1:length(nVec)

        N = nVec(n);

        % ANA-II1
        indicesAug = AugNestedII1Gen(N);
        DOFmaxAug(n) = max(indicesAug);

        % MISC
        indicesMISC = MISC(N);
        DOFmaxMISC(n) = max(indicesMISC);

    end

    % MRA
    DOFmaxOptTemp       = accumarray(OptResult(:, 2), OptResult(:, 1), [], @max);
    DOFmaxOpt           = DOFmaxOptTemp(DOFmaxOptTemp > 0);
    NsOpt               = unique(OptResult(:, 2));

    % SMART
    DOFmaxE_SMARTTemp   = accumarray(SMART_Result(:, 2), SMART_Result(:, 1), [], @max);
    DOFmaxSMART        = DOFmaxE_SMARTTemp(DOFmaxE_SMARTTemp > 0);

    % E-SMART
    DOFmaxE_SMARTTemp  = accumarray(E_SMART_Result(:, 2), E_SMART_Result(:, 1), [], @max);
    DOFmaxE_SMART      = DOFmaxE_SMARTTemp(DOFmaxE_SMARTTemp > 0);

    % CCA
    maxDOFCCA           = arrayfun(@(i) max(CCAResult(CCAResult(:,2) == i, 1)), nVec(1):nMax);

    % - Plot Results
    figure('Name','Number of Sensors','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/3) floor(ScreenSize(4)/2)]);
    movegui('west')
    hold on
    minSensor = min([NsSMART; NsE_SMART(NsE_SMART~=0); NsAug(NsAug~=0); NsMISC(NsMISC~=0); OptResult(find(OptResult(:,1) == 2): find(OptResult(:,1) == DOFmax), 2); CCAResult_prime_free(1:find(CCAResult_prime_free(:,1) == DOFmax),2)]);
    maxSensor = max([NsSMART; NsE_SMART(NsE_SMART~=0); NsAug(NsAug~=0); NsMISC(NsMISC~=0); OptResult(find(OptResult(:,1) == 2): find(OptResult(:,1) == DOFmax), 2); CCAResult_prime_free(1:find(CCAResult_prime_free(:,1) == DOFmax),2)]);
    yStep = 0.7/5;
    for i = minSensor:maxSensor
        X = [0 DOFmax DOFmax 0];
        Y = [i-0.5 i-0.5 i+0.5 i+0.5];

        fill(X,Y,[1,1,1]./(1.5^(mod(i,2))), 'EdgeColor','none','FaceAlpha',0.5)
    end
    h(1) = plot(2:DOFmax, NsSMART-5/2*yStep, '*', 'MarkerSize', 6, 'LineWidth', 1.25);
    h(2) = plot(5:DOFmax, NsE_SMART(NsE_SMART~=0)-3/2*yStep, 'v','MarkerSize', 6, 'LineWidth', 1.25);
    h(3) = plot(DOFaug(DOFaug ~= 0), NsAug(NsAug~=0)-1/2*yStep, 'd','MarkerSize', 6, 'LineWidth', 1.25);
    h(4) = plot(DOFmisc(DOFmisc ~= 0), NsMISC(NsMISC~=0)+1/2*yStep, 'p','MarkerSize', 6, 'LineWidth', 1.25);
    h(5) = plot(CCAResult_prime_free(1:find(CCAResult_prime_free(:,1) == DOFmax),1), CCAResult_prime_free(find(CCAResult_prime_free(:,1) == 4): find(CCAResult_prime_free(:,1) == DOFmax), 2)+3/2*yStep, '+','MarkerSize', 6, 'LineWidth', 1.25);
    h(6) = plot(2:DOFmax, OptResult(find(OptResult(:,1) == 2): find(OptResult(:,1) == DOFmax), 2)+5/2*yStep, 's', 'MarkerSize', 6, 'LineWidth', 1.25);
    xlabel('Array Length ($N_A$)','Interpreter','latex','FontSize',14')
    ylabel('Number of Sensors ($N_s$)','Interpreter','latex','FontSize',14')
    legend(h(1:6),'SMART', 'E-SMART', 'ANA-II1', 'MISC', 'CCA', 'MRA', 'Interpreter','latex','Location','nw','FontSize',12)
    xlim([0 DOFmax])
    xticks(0:4:32)
    yticks(3:1:15)
    ylim([minSensor-0.5 maxSensor+0.5])
    grid on
    ax = gca;
    ax.YGrid = 'off';
    ax.TickDir = 'none';
    ax.FontSize = 12;
    axis square
    box on

    figure('Name','Diff Number of Sensors','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/3) floor(ScreenSize(4)/2)]);
    movegui('center')
    hold on
    OptNs = OptResult(find(OptResult(:,1) == 2): find(OptResult(:,1) == DOFmax), 2);
    minDiff = min([NsSMART-OptNs; NsE_SMART(NsE_SMART ~= 0)-OptNs(4:end); NsAug(NsAug~=0)-OptNs(DOFaug(DOFaug ~= 0)-1); NsMISC(NsMISC~=0)-OptNs(DOFmisc(DOFmisc ~= 0)-1); CCAResult_prime_free(find(CCAResult_prime_free(:,1) == 4): find(CCAResult_prime_free(:,1) == DOFmax), 2) - OptNs(CCAResult_prime_free(1:find(CCAResult_prime_free(:,1) == DOFmax),1) - 2)]);
    maxDiff = max([NsSMART-OptNs; NsE_SMART(NsE_SMART ~= 0)-OptNs(4:end); NsAug(NsAug~=0)-OptNs(DOFaug(DOFaug ~= 0)-1); NsMISC(NsMISC~=0)-OptNs(DOFmisc(DOFmisc ~= 0)-1); CCAResult_prime_free(find(CCAResult_prime_free(:,1) == 4): find(CCAResult_prime_free(:,1) == DOFmax), 2) - OptNs(CCAResult_prime_free(1:find(CCAResult_prime_free(:,1) == DOFmax),1) - 2)]);
    yStep = 1.4/4;
    for i = minDiff:maxDiff
        X = [0 DOFmax DOFmax 0];
        Y = [i-0.5 i-0.5 i+0.5 i+0.5];

        fill(X,Y,[1,1,1]./(1.5^(mod(i,2))), 'EdgeColor','none','FaceAlpha',0.5)
    end
    h(7) = plot(2:DOFmax,NsSMART - OptNs - 1*yStep, '*', 'MarkerSize', 6, 'LineWidth', 1.25);
    h(8) = plot(5:DOFmax,NsE_SMART(NsE_SMART ~= 0) - OptNs(4:end) -1/2*yStep, 'v','MarkerSize', 6, 'LineWidth', 1.25);
    h(9) = plot(DOFaug(DOFaug ~= 0), NsAug(NsAug~=0) - OptNs(DOFaug(DOFaug ~= 0)-1) -0*yStep, 'd','MarkerSize', 6, 'LineWidth', 1.25);
    h(10) = plot(DOFmisc(DOFmisc ~= 0), NsMISC(NsMISC~=0) - OptNs(DOFmisc(DOFmisc ~= 0)-1) +1/2*yStep, 'p','MarkerSize', 6, 'LineWidth', 1.25);
    h(11) = plot(CCAResult_prime_free(1:find(CCAResult_prime_free(:,1) == DOFmax),1), CCAResult_prime_free(find(CCAResult_prime_free(:,1) == 4): find(CCAResult_prime_free(:,1) == DOFmax), 2) - OptNs(CCAResult_prime_free(1:find(CCAResult_prime_free(:,1) == DOFmax),1) - 2)+ 1*yStep, '+','MarkerSize', 6, 'LineWidth', 1.25);
    xlabel('Array Length ($N_A$)','Interpreter','latex','FontSize',12')
    ylabel('Difference','Interpreter','latex','FontSize',12')
    legend(h(7:11),'SMART', 'E-SMART', 'ANA-II1', 'MISC', 'CCA', 'Interpreter','latex','Location','nw','FontSize',12)
    xlim([0 DOFmax])
    xticks(0:4:32)
    ylim([minDiff-0.5 maxDiff+0.5])
    grid on
    ax = gca;
    ax.YGrid = 'off';
    ax.TickDir = 'none';
    ax.FontSize = 12;
    axis square
    box on

    figure('Name','Maximum Array Length','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/3) floor(ScreenSize(4)/2)]);
    movegui('east')
    plot(nVec,DOFmaxSMART(1:length(nVec)), '*', 'MarkerSize', 6, 'LineWidth', 1.25), hold on
    plot(nVec(2:end), DOFmaxE_SMART(1:length(nVec)-1), 'v','MarkerSize', 8, 'LineWidth', 1.25)
    plot(nVec,DOFmaxAug, 'd','MarkerSize', 10, 'LineWidth', 1.25)
    plot(nVec,DOFmaxMISC, 'p','MarkerSize', 10, 'LineWidth', 1.25)
    plot(nVec,maxDOFCCA, '+','MarkerSize', 12,'LineWidth', 1.25)
    plot(NsOpt(2:end), DOFmaxOpt(2:end), 's', 'MarkerSize', 10, 'LineWidth', 1.25)
    legend('SMART', 'E-SMART', 'ANA-II1', 'MISC', 'CCA', 'MRA', 'Interpreter','latex','Location','nw','FontSize',11)
    xlabel('Number of Sensors ($N_s$)','Interpreter','latex','FontSize',12')
    ylabel('Maximum Achieved Array Length','Interpreter','latex','FontSize',12')
    xlim([0 nMax])
    xticks(0:4:nMax)
    grid on
    ax = gca;
    ax.FontSize = 12;
    axis square
    box on

    drawnow
    clearvars -except figureNum ScreenSize

end

%%  Figure 6 (Leakage vs Ns - Approximation Coupling)

if any(6 == figureNum)

    % - Array Concfiguration
    load SMART_Result
    load E_SMART_Result
    load CCAresultTempUnsort
    load MRAconfig

    DOFSMARTTemp   = accumarray(SMART_Result(:, 2), SMART_Result(:, 1), [], @max);
    DOFSMART       = DOFSMARTTemp(DOFSMARTTemp > 0);
    DOFESMARTTemp   = accumarray(E_SMART_Result(:, 2), E_SMART_Result(:, 1), [], @max);
    DOFESMART       = DOFESMARTTemp(DOFESMARTTemp > 0);

    % - Number of Sensors
    nVec            = 15:75;

    % - Coupling Parameter
    B               = inf;

    % - Pre Allocation
    cnt             = 0;

    Lmra            = zeros(1,length(nVec));
    Lsmart         = zeros(1,length(nVec));
    Lesmart        = zeros(1,length(nVec));
    LAug            = zeros(1,length(nVec));
    Lmisc           = zeros(1,length(nVec));
    Lcca            = zeros(1,length(nVec));
    Lula            = zeros(1,length(nVec));

    % - Main Loop

    for N = nVec

        cnt = cnt + 1;

        CanCCA = CCAresultTempUnsort(CCAresultTempUnsort(:,5) == N,:);
        [~, idx] = max(CanCCA(:, 1));
        k = CanCCA(idx,2);
        m = CanCCA(idx,3);
        n = CanCCA(idx,4);

        indicesCCA = CCA(k,m,n);
        indicesULA = 0:N-1;
        if N <= 25
            indicesMRA = MRAconfig{N};
            Cmra = cFunction(abs(indicesMRA - indicesMRA'), B);
            Lmra(cnt) = norm(Cmra - eye(N), 'fro') / norm(Cmra, 'fro');
        end

        indicesSMART = SMART(DOFSMART(N-2),'result','off');
        indicesESMART = E_SMART(DOFESMART(N-3),'result','off');
        indicesMISC = MISC(N);
        indicesAug = AugNestedII1Gen(N);

        Csmart = cFunction(abs(indicesSMART - indicesSMART'), B);
        Cesmart = cFunction(abs(indicesESMART - indicesESMART'), B);
        Caug = cFunction(abs(indicesAug - indicesAug'), B);
        Cmisc = cFunction(abs(indicesMISC - indicesMISC'), B);
        Ccca = cFunction(abs(indicesCCA - indicesCCA'), B);
        CULA = cFunction(abs(indicesULA - indicesULA'), B);

        % - Coupling Leakage
        Lsmart(cnt) = norm(Csmart - eye(N), 'fro') / norm(Csmart, 'fro');
        Lesmart(cnt) = norm(Cesmart - eye(N), 'fro') / norm(Cesmart, 'fro');
        LAug(cnt) = norm(Caug - eye(N), 'fro') / norm(Caug, 'fro');
        Lmisc(cnt) = norm(Cmisc - eye(N), 'fro') / norm(Cmisc, 'fro');
        Lcca(cnt) = norm(Ccca - eye(N), 'fro') / norm(Ccca, 'fro');
        Lula(cnt) = norm(CULA - eye(N), 'fro') / norm(CULA, 'fro');

    end

    % - Plot Results

    figure('Name','Leakave vs Number of Sensors','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/3) floor(ScreenSize(4)/4)]);
    movegui('center')
    semilogy(nVec,Lsmart,'LineWidth',1.5), hold on
    semilogy(nVec,Lesmart,'LineWidth',1.5)
    semilogy(nVec,LAug,'LineWidth',1.5)
    semilogy(nVec,Lmisc,'LineWidth',1.5)
    semilogy(nVec,Lcca,'LineWidth',1.5)
    semilogy(nVec(nVec <= 25), Lmra(1:max([find(nVec <= 25), 0])), '--','LineWidth',1.5)
    semilogy(nVec,Lula,'LineWidth',1.5)
    legend({'SMART', 'E-SMART', 'ANA-II1', 'MISC', 'CCA', 'MRA', 'ULA'}, ...
        'Interpreter','latex', 'FontSize',10, 'Location','northoutside', 'NumColumns', 4);
    xlabel('$N_s$','Interpreter','Latex','FontSize',12)
    ylabel('$L$','Interpreter','Latex','FontSize',12)
    grid minor
    xlim([nVec(1) nVec(end)])
    xticks(15:5:75)
    ylim padded

    drawnow
    clearvars -except figureNum ScreenSize

end

%%  Figure 7 (Leakage and Ns vs DOF)

if any(7 == figureNum)

    % - DOF Range
    DOFvec          = 20:20:1e3;

    % - Coupling Parameter
    B               = inf;

    % - Pre Allocation
    cnt             = 0;
    NSsmart        = zeros(1,length(DOFvec));
    NSesmart        = zeros(1,length(DOFvec));
    Lsmart         = zeros(1,length(DOFvec));
    Lesmart         = zeros(1,length(DOFvec));

    % - Main Loop

    for DOF = DOFvec

        cnt = cnt + 1;

        % Array Configuration
        indicesSMART = SMART(DOF,'result','off');
        indicesESMART = E_SMART(DOF,'result','off');

        Csmart = cFunction(abs(indicesSMART - indicesSMART'), B);
        Cesmart = cFunction(abs(indicesESMART - indicesESMART'), B);

        NSsmart(cnt) = length(indicesSMART);
        NSesmart(cnt) = length(indicesESMART);

        Lsmart(cnt) = norm(Csmart - eye(length(Csmart)), 'fro') / norm(Csmart, 'fro');
        Lesmart(cnt) = norm(Cesmart - eye(length(Cesmart)), 'fro') / norm(Cesmart, 'fro');

    end

    % - Plot Results

    figure('Name','Leakage vs DOF','NumberTitle','off', ...
        'Position', [0 0 floor(ScreenSize(3)/3) floor(ScreenSize(4)/3)]);
    movegui('center')
    yyaxis left
    semilogy(DOFvec, Lsmart, 'LineWidth', 1.5), hold on
    semilogy(DOFvec, Lesmart, 'LineWidth', 1.5)
    ylabel('$L$', 'Interpreter','latex','FontSize',12)
    yyaxis right
    plot(DOFvec, NSsmart, 'LineWidth', 1.5)
    hold on
    plot(DOFvec, NSesmart, 'LineWidth', 1.5)
    ylabel('$N_s$', 'Interpreter','latex','FontSize',12)
    grid minor
    xlabel('DOF', 'Interpreter','latex','FontSize',12)
    xlim([20 1e3])
    xticks([20 100:100:1e3])
    ylim padded
    legend({'SMART (Leakage)', 'E-SMART (Leakage)', 'SMART ($N_s$)', 'E-SMART ($N_s$)'}, ...
        'Interpreter','latex', 'FontSize',10, 'Location','n');

    drawnow
    clearvars -except figureNum ScreenSize

end

%%  Figure 8 (MUSIC Spectra)

if any(8 == figureNum)

    N               = 11;
    DOA             = linspace(-20,20,N);
    nSnapshots      = 500;
    nSig            = length(DOA);
    Nr              = 14;
    SNR             = 0;

    % - Array Configuration
    dVecCCA         = CCA(8,1,7)';
    dVecMRA         = [0, 1, 2, 8, 15, 16, 26, 36, 46, 56, 59, 63, 65, 68]';
    dVecSMART      = SMART(58,'result','off')';
    dVecE_SMART     = E_SMART(55,'result','off')';
    dVecMISC        = MISC(14)';
    dVecAug         = AugNestedII1Gen(14)';

    dMat            = [dVecSMART dVecE_SMART dVecAug dVecMISC dVecCCA dVecMRA];

    % - Main Loop
    % Pre Allocation
    theta           = linspace(-90,90,1e4);
    doas            = cell(length(dMat),1);
    spec            = zeros(length(theta),length(dMat));

    for method = 1:size(dMat,2)

        dVec = dMat(:,method);

        signals = exp(1j * 2 * pi * rand(nSig, nSnapshots));

        arrayResponse = exp(-1j * pi * (0:(max(dVec)-1))' * sind(DOA(:)'));
        receivedSignals = awgn(arrayResponse * signals, SNR);

        covmat = (receivedSignals * receivedSignals') / nSnapshots;

        [doas{method}, spec(:,method)] = musicdoa(covmat, nSig, 'ScanAngles', theta);

    end

    % - Results
    % SMART
    RMSESMART      = sqrt(1/N * sum((DOA - sort(doas{1})).^2));

    figure('Name','MUSIC Spectrum SMART','NumberTitle','off',...
        'Position', [0 2*floor(ScreenSize(4)/5) floor(ScreenSize(3)/6) floor(ScreenSize(4)/9)]);
    xline(DOA, 'LineWidth', 1), hold on
    plot(theta, 20*log10(spec(:,1)/max(spec(:,1))), 'LineWidth', 1.5)
    set(gca, 'YMinorGrid', 'off');
    set(gca, 'YMinorTick', 'off');
    xlabel('$\bar{\theta}$ (degree $^\circ$)','Interpreter','latex','FontSize',10')
    ylabel('$P(\bar{\theta})$ (dB)','Interpreter','latex','FontSize',10')
    xlim([-90 90])
    xticks(-90:30:90)
    box on
    grid on

    % E-SMART
    RMSEE_SMART     = sqrt(1/N * sum((DOA - sort(doas{2})).^2));

    figure('Name','MUSIC Spectrum E-SMART','NumberTitle','off',...
        'Position', [floor(ScreenSize(3)/6) 2*floor(ScreenSize(4)/5) floor(ScreenSize(3)/6) floor(ScreenSize(4)/9)]);
    xline(DOA, 'LineWidth', 1), hold on
    plot(theta, 20*log10(spec(:,2)/max(spec(:,2))), 'LineWidth', 1.5)
    set(gca, 'YMinorGrid', 'off');
    set(gca, 'YMinorTick', 'off');
    xlabel('$\bar{\theta}$ (degree $^\circ$)','Interpreter','latex','FontSize',10')
    ylabel('$P(\bar{\theta})$ (dB)','Interpreter','latex','FontSize',10')
    xlim([-90 90])
    xticks(-90:30:90)
    box on
    grid on

    % ANA-II1
    RMSEaug         = sqrt(1/N * sum((DOA - sort(doas{3})).^2));

    figure('Name','MUSIC Spectrum ANA-II1','NumberTitle','off',...
        'Position', [2*floor(ScreenSize(3)/6) 2*floor(ScreenSize(4)/5) floor(ScreenSize(3)/6) floor(ScreenSize(4)/9)]);
    xline(DOA, 'LineWidth', 1), hold on
    plot(theta, 20*log10(spec(:,3)/max(spec(:,3))), 'LineWidth', 1.5)
    set(gca, 'YMinorGrid', 'off');
    set(gca, 'YMinorTick', 'off');
    xlabel('$\bar{\theta}$ (degree $^\circ$)','Interpreter','latex','FontSize',10')
    ylabel('$P(\bar{\theta})$ (dB)','Interpreter','latex','FontSize',10')
    xlim([-90 90])
    xticks(-90:30:90)
    box on
    grid on

    % MISC
    RMSEmisc        = sqrt(1/N * sum((DOA - sort(doas{4})).^2));

    figure('Name','MUSIC Spectrum MISC','NumberTitle','off',...
        'Position', [3*floor(ScreenSize(3)/6) 2*floor(ScreenSize(4)/5) floor(ScreenSize(3)/6) floor(ScreenSize(4)/9)]);
    xline(DOA, 'LineWidth', 1), hold on
    plot(theta, 20*log10(spec(:,4)/max(spec(:,4))), 'LineWidth', 1.5)
    set(gca, 'YMinorGrid', 'off');
    set(gca, 'YMinorTick', 'off');
    xlabel('$\bar{\theta}$ (degree $^\circ$)','Interpreter','latex','FontSize',10')
    ylabel('$P(\bar{\theta})$ (dB)','Interpreter','latex','FontSize',10')
    xlim([-90 90])
    xticks(-90:30:90)
    box on
    grid on

    % CCA
    RMSEcca         = sqrt(1/N * sum((DOA - sort(doas{5})).^2));

    figure('Name','MUSIC Spectrum CCA','NumberTitle','off',...
        'Position', [4*floor(ScreenSize(3)/6) 2*floor(ScreenSize(4)/5) floor(ScreenSize(3)/6) floor(ScreenSize(4)/9)]);
    xline(DOA, 'LineWidth', 1), hold on
    plot(theta, 20*log10(spec(:,5)/max(spec(:,5))), 'LineWidth', 1.5)
    set(gca, 'YMinorGrid', 'off');
    set(gca, 'YMinorTick', 'off');
    xlabel('$\bar{\theta}$ (degree $^\circ$)','Interpreter','latex','FontSize',10')
    ylabel('$P(\bar{\theta})$ (dB)','Interpreter','latex','FontSize',10')
    xlim([-90 90])
    xticks(-90:30:90)
    box on
    grid on

    % MRA
    RMSEmra         = sqrt(1/N * sum((DOA - sort(doas{6})).^2));

    figure('Name','MUSIC Spectrum MRA','NumberTitle','off',...
        'Position', [5*floor(ScreenSize(3)/6) 2*floor(ScreenSize(4)/5) floor(ScreenSize(3)/6) floor(ScreenSize(4)/9)]);
    xline(DOA, 'LineWidth', 1), hold on
    plot(theta, 20*log10(spec(:,6)/max(spec(:,6))), 'LineWidth', 1.5)
    set(gca, 'YMinorGrid', 'off');
    set(gca, 'YMinorTick', 'off');
    xlabel('$\bar{\theta}$ (degree $^\circ$)','Interpreter','latex','FontSize',10')
    ylabel('$P(\bar{\theta})$ (dB)','Interpreter','latex','FontSize',10')
    xlim([-90 90])
    xticks(-90:30:90)
    box on
    grid on

    fprintf('========================= RMSE =========================\n');
    fprintf('SMART                = %.4f\n', RMSESMART);
    fprintf('E-SMART              = %.4f\n', RMSEE_SMART);
    fprintf('ANA-II1              = %.4f\n', RMSEaug);
    fprintf('MISC                 = %.4f\n', RMSEmisc);
    fprintf('CCA                  = %.4f\n', RMSEcca);
    fprintf('MRA                  = %.4f\n', RMSEmra);
    fprintf('========================================================\n');

    drawnow
    clearvars -except figureNum ScreenSize

end

%%  Figure 9 (Comparison with Non-Hole-Free Array)

if any(9 == figureNum)

    % - (a)

    Ns              = 5;
    indicesCell{1}  = 2.^(0:Ns) - 1;
    indicesCell{2}  = E_SMART(31,'result','off');

    theta0Deg       = 0;
    thetaDeg        = -90:0.25:90;
    k               = 2*pi;

    % Pre Allocation
    psl             = zeros(length(indicesCell),1);
    p               = gobjects(1,2);

    figure('Name','Beampattern Response','NumberTitle','off');
    hold on
    movegui('center')

    % - Main loop

    for i = 1:length(indicesCell)

        indices = indicesCell{i};

        differences = [0; sort(unique(abs(diff(nchoosek(indices, 2), 1, 2))))];

        pos = differences * 0.5;

        w = [];

        if isempty(w)
            w = ones(size(pos));
        end

        theta = deg2rad(thetaDeg);
        theta0 = deg2rad(theta0Deg);

        phase = (sin(theta) - sin(theta0));
        AF = zeros(size(theta));

        for n = 1:numel(pos)
            AF = AF + w(n) .* exp(1j * k * pos(n) .* phase);
        end

        AF = AF / max(abs(AF));
        AFdB = 20*log10(abs(AF));
        psl(i) = sidelobelevel(AFdB);

        p(i) = plot(thetaDeg, AFdB);

    end

    grid minor
    xlabel('$\theta$', 'Interpreter', 'Latex');
    ylabel('Beampattern Response (dB)', 'Interpreter', 'Latex');
    ylim([-40 0])
    xlim([-90 90])
    xticks(-90:30:90)
    legend(p,{'$2$-ELA', 'E-SMART'}, 'Interpreter', 'Latex')
    box on
    yl = yline(psl,'--','HandleVisibility','off');

    for i = 1:numel(p)

        yl(i).Color = p(i).Color;

        str = sprintf('$%.2f\\,\\mathrm{dB}$', yl(i).Value);

        text(-70, yl(i).Value+1.5, str, ...
            'Interpreter','latex', ...
            'Color', p(i).Color, ...
            'FontSize', 12);

    end

    % - (b)

    Ns              = 5;
    d               = 0.5;
    DOA             = sort(-25:5:25);
    KTrue           = length(DOA);
    SNRdB           = 0;
    Nsnap           = 500;
    ang             = -90:0.1:90;

    indicesCell{1}  = E_SMART(31,'result','off');
    indicesCell{2}  = (2.^(0:Ns)-1);

    pkN = @(p,x) findpeaks(p, x, 'SortStr','descend');

    for cnt = 1:length(indicesCell)

        indices = indicesCell{cnt};

        M = numel(indices);

        pos = indices * d;

        D = abs(indices(:) - indices(:).');
        D = unique(D(:)).';

        L = 0;

        while ismember(L+1, D)
            L = L+1;
        end

        k0 = 2*pi;
        A  = exp(1j*k0*(pos.').*sind(DOA));

        S = (randn(KTrue,Nsnap)+1j*randn(KTrue,Nsnap))/sqrt(2);
        S = S * 10^(SNRdB/20);
        N = (randn(M,Nsnap)+1j*randn(M,Nsnap))/sqrt(2);
        X = A*S + N;

        Rxx = (X*X')/Nsnap;
        Rxx = (Rxx + Rxx')/2;

        diffMat = indices(:) - indices(:).';
        kset = -L:L;
        rhat = zeros(1, numel(kset));
        for ii = 1:numel(kset)
            k = kset(ii);
            mask = (diffMat == k);
            vals = Rxx(mask);
            rhat(ii) = mean(vals);
        end

        rhat = (rhat + conj(fliplr(rhat)))/2;

        Rvss = toeplitz(rhat(L+1:end), conj(rhat(L+1:-1:1)));
        Rvss = (Rvss + Rvss')/2;

        Kss = min(KTrue, L);
        [Uss, Dss] = eig(Rvss);
        [~,ord] = sort(real(diag(Dss)),'descend');
        Uss = Uss(:,ord);
        Ens = Uss(:, Kss+1:end);

        posVss = (0:L)*d;
        aVss = @(th) exp(1j*k0*(posVss.').*sind(th));

        Pss = zeros(size(ang));
        for i = 1:numel(ang)

            av = aVss(ang(i));
            Pss(i) = 1./real(av'*(Ens*Ens')*av);

        end

        Pss = Pss./max(Pss);
        PssdB = 10*log10(Pss);

        figure('Name','MUSIC Spectrum','NumberTitle','off',...
            'Position', [(cnt-1)*(ScreenSize(3) - floor(ScreenSize(3)/4)) 2*floor(ScreenSize(4)/5) floor(ScreenSize(3)/4) floor(ScreenSize(4)/5)]);
        xline(DOA, 'LineWidth', 1), hold on
        plot(ang, PssdB, 'LineWidth',1.2); grid on
        xlabel('$\bar{\theta}$', 'Interpreter', 'Latex');
        ylabel('$P(\bar{\theta})$ (dB)', 'Interpreter', 'Latex');
        xlim([-90 90])
        xticks(-90:30:90)
        box on

    end

    drawnow
    clearvars -except figureNum ScreenSize

end

%%  Figure 10 (Probability of Resolution)

if any(10 == figureNum)

    % - Taget and MUSIC Parameters
    nSnapShot       = 1e3;
    thetaScan       = -90:0.1:90;
    tSpacing        = 5;
    theta           = [-tSpacing/2 tSpacing/2];
    K               = length(theta);

    aSpacing        = 0.5;

    % - Coupling Parameter
    B               = inf;

    % - Simulation Parameters
    snrVec          = -10:1:25;
    nMC             = 1e3;

    % - Array configuration
    load MRAconfig
    N               = 20;
    indicesCell     = {
        SMART(112, 'result', 'off');
        E_SMART(105, 'result', 'off');
        AugNestedII1Gen(N);
        MISC(N);
        MRAconfig{N};
        };
    methodName      = {'SMART'
        'E-SMART'
        'ANA-II1'
        'MISC'
        'MRA'};

    % - Pre Allocation
    Pres            = zeros(1, length(snrVec));

    % - Parfor Pre Initialization
    poolobj         = gcp('nocreate');
    if isempty(poolobj)
        c = parcluster;
        c.NumWorkers = 8;
        parpool(c)
    end

    % - Main Loop

    for method = 1:length(indicesCell)

        fprintf('Method: %s ...\n', methodName{method});

        indices = indicesCell{method};

        C = cFunction(abs(indices - indices'), B);

        A = exp(1j * 2 * pi * aSpacing * (indices).' * sind(theta));

        for iSNR = 1:length(snrVec)

            snr = snrVec(iSNR);
            pr = 0;

            parfor iMC = 1:nMC
                signals = exp(1j * 2 * pi * rand(K, nSnapShot));
                x = (C * A * signals).';
                RxxFull = (x' * x) / nSnapShot;

                PosIdxReconstruct = -max(indices):max(indices);
                PosIdxDiff = indices-indices.';
                [~, CovIdx] = ismember(PosIdxReconstruct(:),PosIdxDiff);

                X = RxxFull(CovIdx);
                Y = awgn(X, snr - 10*log10(2*length(indices)-1),'measured');

                covmat = (Y * Y');

                ssfactor = (length(X)+1)/2;
                Rss = spsmooth(covmat,ssfactor);

                thetaHatTemp = musicdoa(Rss, K, 'ScanAngles', thetaScan);
                thetaHat = sort(thetaHatTemp(:).');

                if all(abs(thetaHat - theta) <= tSpacing/5)
                    pr(iMC) = 1;
                end

            end

            Pres(method,iSNR) = sum(pr) / nMC;

        end

    end

    % - Plot Results

    figure('Name','Probability of Resolution','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/3) floor(ScreenSize(3)/3)]);
    movegui('center')
    plot(snrVec, Pres, 'LineWidth', 2)
    grid on
    xlabel('SNR (dB)', 'Interpreter', 'latex','FontSize', 12)
    ylabel('Probability of Resolution', 'Interpreter', 'latex','FontSize', 12)
    legend('SMART','E-SMART','ANA-II1','MISC','MRA', 'Interpreter', 'latex','FontSize', 12, 'Location', 'e')

    drawnow
    clearvars -except figureNum ScreenSize

end

%%  Figure 11 (RMSE Monte-Carlo)

if any(11 == figureNum)

    % Parameters
    Ns                  = 14;

    nSig                = 20;
    DOA                 = linspace(-72,72,nSig);
    nSnapshots          = 1e3;
    snrVec              = -10:5:30;
    nMonteCarlo         = 1e3;
    theta               = linspace(-90, 90, 1e5);

    % - Array Configuration
    DOFSMART           = 58;
    DOFE_SMART          = 59;

    indicesSMART       = SMART(DOFSMART,'result','off');
    indicesE_SMART      = E_SMART(DOFE_SMART,'result','off');
    indicesAug          = AugNestedII1Gen(Ns);
    indicesMISC         = MISC(Ns);
    indicesCCA          = CCA(8,1,7);
    indicesMRA          = [0, 1, 2, 8, 15, 16, 26, 36, 46, 56, 59, 63, 65, 68];

    maxDOFvec           = [max(indicesSMART) max(indicesE_SMART) max(indicesAug) max(indicesMISC) max(indicesCCA) max(indicesMRA)];

    % - Main Loop
    % Pre Allocation
    RMSE                = zeros(nMonteCarlo,1);
    meanRMSE            = zeros(length(maxDOFvec),length(snrVec));

    for method = 1:length(maxDOFvec)

        maxDOF =  maxDOFvec(method);

        for i = 1:length(snrVec)

            SNR = snrVec(i);

            for j = 1:nMonteCarlo

                signals = exp(1j * 2 * pi * rand(nSig, nSnapshots));

                arrayResponse = exp(-1j * pi * (0:(maxDOF))' * sind(DOA));
                receivedSignals = awgn(arrayResponse * signals, SNR);

                covmat = (receivedSignals * receivedSignals') / nSnapshots;

                [doas, spec] = musicdoa(covmat, nSig, 'ScanAngles', theta);

                RMSE(j) = sqrt(1/length(DOA) * sum((DOA - sort(doas)).^2));

            end

            meanRMSE(method,i) = mean(RMSE);

        end

    end

    % - Plot Results
    lineStyles          = {'-*', '-v', '-d', '-p', '-+', '-s'};

    fhdl(1)             = figure('Name','RMSE vs SNR','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/3) floor(ScreenSize(4)/2.5)]);
    movegui('center')
    h1                  = subplot(1,4,1:3);
    hold on
    for i = 1:size(meanRMSE,1)
        semilogy(snrVec,meanRMSE(i,:),lineStyles{i},'LineWidth',2,'MarkerSize',8)
    end
    legend('SMART', 'E-SMART', 'ANA-II1', 'MISC', 'CCA', 'MRA', 'Interpreter','latex','Location','best','FontSize',12)
    ylabel('RMSE (degree)','Interpreter','latex','FontSize',12)
    grid minor
    box on
    rectPos             = [4, 0.0015, 2, 0.004];
    rectangle('Position', ...
        rectPos, ...
        'EdgeColor', '#FF00FF', 'LineWidth', 1.5);

    h2                  = subplot(1,4,4);
    hold on
    for i = 1:size(meanRMSE,1)
        semilogy(snrVec,meanRMSE(i,:),lineStyles{i},'LineWidth',2,'MarkerSize',8)
    end
    xlim([4 6])
    grid minor
    box on

    han                 = axes(fhdl(1),'visible','off');
    han.XLabel.Visible  = 'on';
    xlabel(han,'SNR (dB)','Interpreter','latex','FontSize',12)
    subplotPos1         = h1.Position;
    subplotPos2         = h2.Position;
    arrowStartX         = subplotPos1(1) + subplotPos1(3) / (snrVec(end) - snrVec(1)) * 14;
    arrowStartY1        = subplotPos1(2) + subplotPos1(4) / 0.03 * rectPos(2);
    arrowStartY2        = subplotPos1(2) + subplotPos1(4) / 0.03 * (rectPos(2) + rectPos(4));
    arrowEndX           = subplotPos2(1);
    arrowEndY1          = subplotPos2(2);
    arrowEndY2          = subplotPos2(4) + subplotPos2(2);
    annotation('arrow', [arrowStartX, arrowEndX], [arrowStartY1, arrowEndY1], 'Color', '#FF00FF', 'LineWidth', 1.5);
    annotation('arrow', [arrowStartX, arrowEndX], [arrowStartY2, arrowEndY2], 'Color', '#FF00FF', 'LineWidth', 1.5);

    drawnow
    clearvars -except figureNum ScreenSize

end

%% Figure 12 (VA Design - Case 1)

if any(12 == figureNum)

    load desiredVA1
    incFlag         = false;
    reverseStr      = '';

    % - Configuration
    desiredVA = desiredVA1;
    txPos = SMART(50,'result','off');

    % VA Truncation
    sumH            = sum(desiredVA,2);
    sumV            = sum(desiredVA,1);
    nZHidx          = find(sumH);
    nZVidx          = find(sumV);
    zeroT           = nZHidx(1) - 1;
    zeroB           = length(sumH) - nZHidx(end);
    zeroL           = nZVidx(1) - 1;
    zeroR           = length(sumV) - nZVidx(end);

    desiredVA       = desiredVA(zeroT+1:end-zeroB,zeroL+1:end-zeroR);
    desiredVAvec    = reshape(fliplr(desiredVA.'),[],1);

    % - Antenna Configration
    tx(txPos + 1)   = 1;
    rxSize          = size(desiredVA) - size(tx) + 1;
    minRx           = ceil(nnz(desiredVA) / length(txPos));

    % - Optimization
    nVar            = prod(rxSize);
    coeffMat        = double_blocked(tx,rxSize);
    epsilon         = 1e-2;
    maxIter         = 1000;
    threshold       = 1e-6;

    % Pre Allocation
    rxStarMat       = zeros(nVar,maxIter);
    binConErr       = zeros(1,maxIter);
    objfunValP      = zeros(1,maxIter);

    % Main Loop
    fprintf('================================================================================\n');
    fprintf('Optimizing ...\n');
    for iter = 1:maxIter

        if iter == 1
            C = eye(nVar);
        else
            C = diag(1./(rx + epsilon));
            minRx = minRx + 1;
        end

        cvx_begin quiet

        % cvx_solver gurobi

        variable rx(nVar) nonnegative

        minimize sum(pow_abs(coeffMat*rx - desiredVAvec,2))

        subject to

        ones(1,nVar) * C*rx <= minRx                            %#ok
        rx <= 1 + 1e-3;                                         %#ok
        C * rx <= 1 + 1e-3;                                     %#ok

        cvx_end

        rxStarMat(:,iter) = C * rx;
        binConErr(iter) = 1/nVar * norm(rxStarMat(:,iter) - round(rxStarMat(:,iter)),1);

        rxOptIter = fliplr(reshape(round(rxStarMat(:,iter)),rxSize(2),rxSize(1))).';
        designedVAiter = conv2(tx,rxOptIter);
        designedVAiter(designedVAiter >= 1) = 1;

        objfunValP(iter) = norm(designedVAiter - desiredVA, 'fro');

        msg = sprintf('Iteration %2d: Objective Function Value = %.3f, Binary Con. Error = %2d\n', ...
            iter, objfunValP(iter), binConErr(iter));
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));

        % Stopping Criterion
        if iter >= 15 && max(objfunValP(iter - 14:iter)) - mean(objfunValP(iter - 14:iter)) <= threshold
            break
        end

        if iter >= 20 && objfunValP(iter - 1) < objfunValP(iter)
            incFlag = true;
            break
        end

        if norm(designedVAiter - desiredVA, 'fro') < 1
            break
        end

    end
    fprintf('================================================================================\n');
    fprintf('Optimization Completed\n');

    rxStar          = (fliplr(reshape(round(rxStarMat(:,iter)),rxSize(2),rxSize(1))).');
    if incFlag
        rxStar = (fliplr(reshape(round(rxStarMat(:,iter-1)),rxSize(2),rxSize(1))).');
    end

    % - Results
    [yVA, xVA]      = ind2sub(size(desiredVA),find(desiredVA == 1));

    DesignedVA      = conv2(tx,rxStar);
    DesignedVA(DesignedVA >=1) = 1;
    [ydVA, xdVA]    = ind2sub(size(DesignedVA),find(DesignedVA == 1));

    [yTx, xTx]      = ind2sub(size(tx),find(tx >= 1));
    [yRx, xRx]      = ind2sub(size(rxStar),find(rxStar == 1));

    % Number of Rx
    NrxStar         = sum(nnz(rxStar));

    % Display Results
    fprintf('--------------------------------------------------------------------------------\n');
    fprintf('Number of Tx Elements     = %3d\n', sum(nnz(tx)));
    fprintf('Number of Rx Elements     = %3d\n', NrxStar);
    fprintf('VA Placement Failure Rate = %.1f Percent \n', nnz(DesignedVA - desiredVA) / nnz(desiredVA) * 100);
    fprintf('================================================================================\n');

    % Plot Results
    figure('Name','Tx-Rx Configuration','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/4) floor(ScreenSize(4)/2.5)]);
    movegui('west')
    plot(xTx,yTx,'o','MarkerSize',5,'LineWidth',1.5), hold on
    plot(xRx,yRx + 2,'v','MarkerSize',5,'LineWidth',1.5)
    xlabel('Horizontal (Half Wavelength)','Interpreter','latex','FontSize',12')
    ylabel('Vertical (Half Wavelength)', 'Interpreter', 'latex','FontSize',12')
    legend('Tx Chain', 'Rx Chain', 'Interpreter','latex','Location','northoutside','FontSize',11)
    axis square
    grid on
    xlim([0 max(max(xTx), max(xRx))+1])
    ylim([0 max(max(yTx), max(yRx))+3])

    figure('Name','Virtual Array Configuration','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/4) floor(ScreenSize(4)/2.5)]);
    movegui('east')
    plot(xVA,yVA,'s','MarkerSize',10,'LineWidth',1.5), hold on
    plot(xdVA,ydVA,'.','MarkerSize',12)
    xlabel('Horizontal (Half Wavelength)','Interpreter','latex','FontSize',12')
    ylabel('Vertical (Half Wavelength)', 'Interpreter', 'latex','FontSize',12')
    legend('Desired VA', 'Designed VA', 'Interpreter','latex','Location','northoutside','FontSize',11)
    axis square
    grid on
    xlim([0 max(max(xVA), max(xdVA))+1])
    ylim([0 max(max(yVA), max(ydVA))+1])

    drawnow
    clearvars -except figureNum ScreenSize

end

%% Figure 13 (VA Design - Case 2)

if any(13 == figureNum)

    load desiredVA6
    incFlag         = false;
    reverseStr      = '';

    % - Configuration
    desiredVA = desiredVA6;
    txPos = SMART(10,'result','off');

    % VA Truncation
    sumH            = sum(desiredVA,2);
    sumV            = sum(desiredVA,1);
    nZHidx          = find(sumH);
    nZVidx          = find(sumV);
    zeroT           = nZHidx(1) - 1;
    zeroB           = length(sumH) - nZHidx(end);
    zeroL           = nZVidx(1) - 1;
    zeroR           = length(sumV) - nZVidx(end);

    desiredVA       = desiredVA(zeroT+1:end-zeroB,zeroL+1:end-zeroR);
    desiredVAvec    = reshape(fliplr(desiredVA.'),[],1);

    % - Antenna Configration
    tx(txPos + 1)   = 1;
    rxSize          = size(desiredVA) - size(tx) + 1;
    minRx           = ceil(nnz(desiredVA) / length(txPos));

    % - Optimization
    nVar            = prod(rxSize);
    coeffMat        = double_blocked(tx,rxSize);
    epsilon         = 1e-2;
    maxIter         = 1000;
    threshold       = 1e-6;

    % Pre Allocation
    rxStarMat       = zeros(nVar,maxIter);
    binConErr       = zeros(1,maxIter);
    objfunValP      = zeros(1,maxIter);

    % Main Loop
    fprintf('================================================================================\n');
    fprintf('Optimizing ...\n');
    for iter = 1:maxIter

        if iter == 1
            C = eye(nVar);
        else
            C = diag(1./(rx + epsilon));
            minRx = minRx + 1;
        end

        cvx_begin quiet

        % cvx_solver gurobi

        variable rx(nVar) nonnegative

        minimize sum(pow_abs(coeffMat*rx - desiredVAvec,2))

        subject to

        ones(1,nVar) * C*rx <= minRx                            %#ok
        rx <= 1 + 1e-3;                                         %#ok
        C * rx <= 1 + 1e-3;                                     %#ok

        cvx_end

        rxStarMat(:,iter) = C * rx;
        binConErr(iter) = 1/nVar * norm(rxStarMat(:,iter) - round(rxStarMat(:,iter)),1);

        rxOptIter = fliplr(reshape(round(rxStarMat(:,iter)),rxSize(2),rxSize(1))).';
        designedVAiter = conv2(tx,rxOptIter);
        designedVAiter(designedVAiter >= 1) = 1;

        objfunValP(iter) = norm(designedVAiter - desiredVA, 'fro');

        msg = sprintf('Iteration %2d: Objective Function Value = %.3f, Binary Con. Error = %2d\n', ...
            iter, objfunValP(iter), binConErr(iter));
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));

        % Stopping Criterion
        if iter >= 15 && max(objfunValP(iter - 14:iter)) - mean(objfunValP(iter - 14:iter)) <= threshold
            break
        end

        if iter >= 20 && objfunValP(iter - 1) < objfunValP(iter)
            incFlag = true;
            break
        end

        if norm(designedVAiter - desiredVA, 'fro') < 1
            break
        end

    end
    fprintf('================================================================================\n');
    fprintf('Optimization Completed\n');

    rxStar          = (fliplr(reshape(round(rxStarMat(:,iter)),rxSize(2),rxSize(1))).');
    if incFlag
        rxStar = (fliplr(reshape(round(rxStarMat(:,iter-1)),rxSize(2),rxSize(1))).');
    end

    % - Results
    [yVA, xVA]      = ind2sub(size(desiredVA),find(desiredVA == 1));

    DesignedVA      = conv2(tx,rxStar);
    DesignedVA(DesignedVA >=1) = 1;
    [ydVA, xdVA]    = ind2sub(size(DesignedVA),find(DesignedVA == 1));

    [yTx, xTx]      = ind2sub(size(tx),find(tx >= 1));
    [yRx, xRx]      = ind2sub(size(rxStar),find(rxStar == 1));

    % Number of Rx
    NrxStar         = sum(nnz(rxStar));

    % Display Results
    fprintf('--------------------------------------------------------------------------------\n');
    fprintf('Number of Tx Elements     = %3d\n', sum(nnz(tx)));
    fprintf('Number of Rx Elements     = %3d\n', NrxStar);
    fprintf('VA Placement Failure Rate = %.1f Percent \n', nnz(DesignedVA - desiredVA) / nnz(desiredVA) * 100);
    fprintf('================================================================================\n');

    % Plot Results
    figure('Name','Tx-Rx Configuration','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/4) floor(ScreenSize(4)/2.5)]);
    movegui('west')
    plot(xTx,yTx,'o','MarkerSize',5,'LineWidth',1.5), hold on
    plot(xRx,yRx + 2,'v','MarkerSize',5,'LineWidth',1.5)
    xlabel('Horizontal (Half Wavelength)','Interpreter','latex','FontSize',12')
    ylabel('Vertical (Half Wavelength)', 'Interpreter', 'latex','FontSize',12')
    legend('Tx Chain', 'Rx Chain', 'Interpreter','latex','Location','northoutside','FontSize',11)
    axis square
    grid on
    xlim([0 max(max(xTx), max(xRx))+1])
    ylim([0 max(max(yTx), max(yRx))+3])

    figure('Name','Virtual Array Configuration','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/4) floor(ScreenSize(4)/2.5)]);
    movegui('east')
    plot(xVA,yVA,'s','MarkerSize',10,'LineWidth',1.5), hold on
    plot(xdVA,ydVA,'.','MarkerSize',12)
    xlabel('Horizontal (Half Wavelength)','Interpreter','latex','FontSize',12')
    ylabel('Vertical (Half Wavelength)', 'Interpreter', 'latex','FontSize',12')
    legend('Desired VA', 'Designed VA', 'Interpreter','latex','Location','northoutside','FontSize',11)
    axis square
    grid on
    xlim([0 max(max(xVA), max(xdVA))+1])
    ylim([0 max(max(yVA), max(ydVA))+1])

    drawnow
    clearvars -except figureNum ScreenSize

end

%% Figure 14 (VA Design - Case 3)

if any(14 == figureNum)

    load desiredVA5
    incFlag         = false;
    reverseStr      = '';

    % - Configuration
    desiredVA = desiredVA5;
    txPos = SMART(10,'result','off');

    % VA Truncation
    sumH            = sum(desiredVA,2);
    sumV            = sum(desiredVA,1);
    nZHidx          = find(sumH);
    nZVidx          = find(sumV);
    zeroT           = nZHidx(1) - 1;
    zeroB           = length(sumH) - nZHidx(end);
    zeroL           = nZVidx(1) - 1;
    zeroR           = length(sumV) - nZVidx(end);

    desiredVA       = desiredVA(zeroT+1:end-zeroB,zeroL+1:end-zeroR);
    desiredVAvec    = reshape(fliplr(desiredVA.'),[],1);

    % - Antenna Configration
    tx(txPos + 1)   = 1;
    rxSize          = size(desiredVA) - size(tx) + 1;
    minRx           = ceil(nnz(desiredVA) / length(txPos));

    % - Optimization
    nVar            = prod(rxSize);
    coeffMat        = double_blocked(tx,rxSize);
    epsilon         = 1e-2;
    maxIter         = 1000;
    threshold       = 1e-6;

    % Pre Allocation
    rxStarMat       = zeros(nVar,maxIter);
    binConErr       = zeros(1,maxIter);
    objfunValP      = zeros(1,maxIter);

    % Main Loop
    fprintf('================================================================================\n');
    fprintf('Optimizing ...\n');
    for iter = 1:maxIter

        if iter == 1
            C = eye(nVar);
        else
            C = diag(1./(rx + epsilon));
            minRx = minRx + 1;
        end

        cvx_begin quiet

        % cvx_solver gurobi

        variable rx(nVar) nonnegative

        minimize sum(pow_abs(coeffMat*rx - desiredVAvec,2))

        subject to

        ones(1,nVar) * C*rx <= minRx                            %#ok
        rx <= 1 + 1e-3;                                         %#ok
        C * rx <= 1 + 1e-3;                                     %#ok

        cvx_end

        rxStarMat(:,iter) = C * rx;
        binConErr(iter) = 1/nVar * norm(rxStarMat(:,iter) - round(rxStarMat(:,iter)),1);

        rxOptIter = fliplr(reshape(round(rxStarMat(:,iter)),rxSize(2),rxSize(1))).';
        designedVAiter = conv2(tx,rxOptIter);
        designedVAiter(designedVAiter >= 1) = 1;

        objfunValP(iter) = norm(designedVAiter - desiredVA, 'fro');

        msg = sprintf('Iteration %2d: Objective Function Value = %.3f, Binary Con. Error = %2d\n', ...
            iter, objfunValP(iter), binConErr(iter));
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));

        % Stopping Criterion
        if iter >= 15 && max(objfunValP(iter - 14:iter)) - mean(objfunValP(iter - 14:iter)) <= threshold
            break
        end

        if iter >= 20 && objfunValP(iter - 1) < objfunValP(iter)
            incFlag = true;
            break
        end

        if norm(designedVAiter - desiredVA, 'fro') < 1
            break
        end

    end
    fprintf('================================================================================\n');
    fprintf('Optimization Completed\n');

    rxStar          = (fliplr(reshape(round(rxStarMat(:,iter)),rxSize(2),rxSize(1))).');
    if incFlag
        rxStar = (fliplr(reshape(round(rxStarMat(:,iter-1)),rxSize(2),rxSize(1))).');
    end

    % - Results
    [yVA, xVA]      = ind2sub(size(desiredVA),find(desiredVA == 1));

    DesignedVA      = conv2(tx,rxStar);
    DesignedVA(DesignedVA >=1) = 1;
    [ydVA, xdVA]    = ind2sub(size(DesignedVA),find(DesignedVA == 1));

    [yTx, xTx]      = ind2sub(size(tx),find(tx >= 1));
    [yRx, xRx]      = ind2sub(size(rxStar),find(rxStar == 1));

    % Number of Rx
    NrxStar         = sum(nnz(rxStar));

    % Display Results
    fprintf('--------------------------------------------------------------------------------\n');
    fprintf('Number of Tx Elements     = %3d\n', sum(nnz(tx)));
    fprintf('Number of Rx Elements     = %3d\n', NrxStar);
    fprintf('VA Placement Failure Rate = %.1f Percent \n', nnz(DesignedVA - desiredVA) / nnz(desiredVA) * 100);
    fprintf('================================================================================\n');

    % Plot Results
    figure('Name','Tx-Rx Configuration','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/4) floor(ScreenSize(4)/2.5)]);
    movegui('west')
    plot(xTx,yTx,'o','MarkerSize',5,'LineWidth',1.5), hold on
    plot(xRx,yRx + 2,'v','MarkerSize',5,'LineWidth',1.5)
    xlabel('Horizontal (Half Wavelength)','Interpreter','latex','FontSize',12')
    ylabel('Vertical (Half Wavelength)', 'Interpreter', 'latex','FontSize',12')
    legend('Tx Chain', 'Rx Chain', 'Interpreter','latex','Location','northoutside','FontSize',11)
    axis square
    grid on
    xlim([0 max(max(xTx), max(xRx))+1])
    ylim([0 max(max(yTx), max(yRx))+3])

    figure('Name','Virtual Array Configuration','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/4) floor(ScreenSize(4)/2.5)]);
    movegui('east')
    plot(xVA,yVA,'s','MarkerSize',10,'LineWidth',1.5), hold on
    plot(xdVA,ydVA,'.','MarkerSize',12)
    xlabel('Horizontal (Half Wavelength)','Interpreter','latex','FontSize',12')
    ylabel('Vertical (Half Wavelength)', 'Interpreter', 'latex','FontSize',12')
    legend('Desired VA', 'Designed VA', 'Interpreter','latex','Location','northoutside','FontSize',11)
    axis square
    grid on
    xlim([0 max(max(xVA), max(xdVA))+1])
    ylim([0 max(max(yVA), max(ydVA))+1])

    drawnow
    clearvars -except figureNum ScreenSize

end

%% Figure 15 (VA Design - Case 4)

if any(15 == figureNum)

    incFlag         = false;
    reverseStr      = '';

    % - Configuration
    desiredVA = ones(20,20);
    txPos = SMART(8,'result','off');

    % VA Truncation
    sumH            = sum(desiredVA,2);
    sumV            = sum(desiredVA,1);
    nZHidx          = find(sumH);
    nZVidx          = find(sumV);
    zeroT           = nZHidx(1) - 1;
    zeroB           = length(sumH) - nZHidx(end);
    zeroL           = nZVidx(1) - 1;
    zeroR           = length(sumV) - nZVidx(end);

    desiredVA       = desiredVA(zeroT+1:end-zeroB,zeroL+1:end-zeroR);
    desiredVAvec    = reshape(fliplr(desiredVA.'),[],1);

    % - Antenna Configration
    tx(txPos + 1)   = 1;
    rxSize          = size(desiredVA) - size(tx) + 1;
    minRx           = ceil(nnz(desiredVA) / length(txPos));

    % - Optimization
    nVar            = prod(rxSize);
    coeffMat        = double_blocked(tx,rxSize);
    epsilon         = 1e-2;
    maxIter         = 1000;
    threshold       = 1e-6;

    % Pre Allocation
    rxStarMat       = zeros(nVar,maxIter);
    binConErr       = zeros(1,maxIter);
    objfunValP      = zeros(1,maxIter);

    % Main Loop
    fprintf('================================================================================\n');
    fprintf('Optimizing ...\n');
    for iter = 1:maxIter

        if iter == 1
            C = eye(nVar);
        else
            C = diag(1./(rx + epsilon));
            minRx = minRx + 1;
        end

        cvx_begin quiet

        % cvx_solver gurobi

        variable rx(nVar) nonnegative

        minimize sum(pow_abs(coeffMat*rx - desiredVAvec,2))

        subject to

        ones(1,nVar) * C*rx <= minRx                            %#ok
        rx <= 1 + 1e-3;                                         %#ok
        C * rx <= 1 + 1e-3;                                     %#ok

        cvx_end

        rxStarMat(:,iter) = C * rx;
        binConErr(iter) = 1/nVar * norm(rxStarMat(:,iter) - round(rxStarMat(:,iter)),1);

        rxOptIter = fliplr(reshape(round(rxStarMat(:,iter)),rxSize(2),rxSize(1))).';
        designedVAiter = conv2(tx,rxOptIter);
        designedVAiter(designedVAiter >= 1) = 1;

        objfunValP(iter) = norm(designedVAiter - desiredVA, 'fro');

        msg = sprintf('Iteration %2d: Objective Function Value = %.3f, Binary Con. Error = %2d\n', ...
            iter, objfunValP(iter), binConErr(iter));
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));

        % Stopping Criterion
        if iter >= 15 && max(objfunValP(iter - 14:iter)) - mean(objfunValP(iter - 14:iter)) <= threshold
            break
        end

        if iter >= 20 && objfunValP(iter - 1) < objfunValP(iter)
            incFlag = true;
            break
        end

        if norm(designedVAiter - desiredVA, 'fro') < 1
            break
        end

    end
    fprintf('================================================================================\n');
    fprintf('Optimization Completed\n');

    rxStar          = (fliplr(reshape(round(rxStarMat(:,iter)),rxSize(2),rxSize(1))).');
    if incFlag
        rxStar = (fliplr(reshape(round(rxStarMat(:,iter-1)),rxSize(2),rxSize(1))).');
    end

    % - Results
    [yVA, xVA]      = ind2sub(size(desiredVA),find(desiredVA == 1));

    DesignedVA      = conv2(tx,rxStar);
    DesignedVA(DesignedVA >=1) = 1;
    [ydVA, xdVA]    = ind2sub(size(DesignedVA),find(DesignedVA == 1));

    [yTx, xTx]      = ind2sub(size(tx),find(tx >= 1));
    [yRx, xRx]      = ind2sub(size(rxStar),find(rxStar == 1));

    % Number of Rx
    NrxStar         = sum(nnz(rxStar));

    % Display Results
    fprintf('--------------------------------------------------------------------------------\n');
    fprintf('Number of Tx Elements     = %3d\n', sum(nnz(tx)));
    fprintf('Number of Rx Elements     = %3d\n', NrxStar);
    fprintf('VA Placement Failure Rate = %.1f Percent \n', nnz(DesignedVA - desiredVA) / nnz(desiredVA) * 100);
    fprintf('================================================================================\n');

    % Plot Results
    figure('Name','Tx-Rx Configuration','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/4) floor(ScreenSize(4)/2.5)]);
    movegui('west')
    plot(xTx,yTx,'o','MarkerSize',5,'LineWidth',1.5), hold on
    plot(xRx,yRx + 2,'v','MarkerSize',5,'LineWidth',1.5)
    xlabel('Horizontal (Half Wavelength)','Interpreter','latex','FontSize',12')
    ylabel('Vertical (Half Wavelength)', 'Interpreter', 'latex','FontSize',12')
    legend('Tx Chain', 'Rx Chain', 'Interpreter','latex','Location','northoutside','FontSize',11)
    axis square
    grid on
    xlim([0 max(max(xTx), max(xRx))+1])
    ylim([0 max(max(yTx), max(yRx))+3])

    figure('Name','Virtual Array Configuration','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/4) floor(ScreenSize(4)/2.5)]);
    movegui('east')
    plot(xVA,yVA,'s','MarkerSize',10,'LineWidth',1.5), hold on
    plot(xdVA,ydVA,'.','MarkerSize',12)
    xlabel('Horizontal (Half Wavelength)','Interpreter','latex','FontSize',12')
    ylabel('Vertical (Half Wavelength)', 'Interpreter', 'latex','FontSize',12')
    legend('Desired VA', 'Designed VA', 'Interpreter','latex','Location','northoutside','FontSize',11)
    axis square
    grid on
    xlim([0 max(max(xVA), max(xdVA))+1])
    ylim([0 max(max(yVA), max(ydVA))+1])

    drawnow
    clearvars -except figureNum ScreenSize

end

%%  Figure 16 (Filled VA Design - 50)

if any(16 == figureNum)

    reverseStr      = '';
    lengthVA        = 50;
    txPos           = E_SMART(25,'result','off');

    % - Antenna Configration
    tx(txPos + 1)   = 1;
    lengthRx        = lengthVA - length(tx) + 1;

    desiredVA       = ones(lengthVA,1);
    minRx           = ceil(lengthVA / length(txPos));

    % - Optimization
    nVar            = lengthRx;
    coeffMat        = double_blocked(tx,[1 lengthRx]);
    epsilon         = 1e-2;
    threshold       = 1e-6;

    subIter         = 0;
    maxIter         = 1e3;
    subIterMax      = 1e2;

    % Pre Allocation
    objfunValP      = zeros(subIterMax,1);
    rxStarMat       = zeros(nVar,maxIter);

    % Main Loop
    fprintf('================================================================================\n');
    fprintf('Optimizing ...\n');
    for iter = 1:maxIter

        subIter = subIter + 1;

        if subIter == 1
            C = eye(nVar);
        else
            C = diag(1./(rx + epsilon));
        end

        cvx_begin quiet

        cvx_solver gurobi

        variable rx(nVar) nonnegative

        minimize sum(pow_abs(coeffMat*rx - desiredVA,2))

        subject to

        ones(1,nVar) * C*rx <= minRx                            %#ok
        rx <= 1 + 1e-3;                                         %#ok
        C * rx <= 1 + 1e-3;                                     %#ok

        cvx_end

        rxStarMat(:,subIter) = C * rx;

        rxOptIter = round(rxStarMat(:,subIter)).';
        designedVAiter = conv2(tx,rxOptIter);
        designedVAiter(designedVAiter >= 1) = 1;

        objfunValP(subIter) = norm(designedVAiter - desiredVA.', 1);

        msg = sprintf('Sub Iteration %2d: Objective Function Value = %3d, Minimum Rx Sensors = %2d \n', ...
            subIter, objfunValP(subIter), minRx);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));

        if norm(designedVAiter - desiredVA.', 1) < 1
            break
        end

        if subIter >= 5 && abs(max(objfunValP(subIter:-1:subIter-4)) - mean(objfunValP(subIter-4:subIter))) <= threshold
            minRx = minRx + 1;
            subIter = 0;
        end

    end
    fprintf('================================================================================\n');
    fprintf('Optimization Completed\n');

    rxStar          = round(rxStarMat(:,subIter)).';

    % - Results
    [yVA, xVA]      = ind2sub(size(desiredVA.'),find(desiredVA.' == 1));

    DesignedVA      = conv2(tx,rxStar);
    DesignedVA(DesignedVA >=1) = 1;
    [ydVA, xdVA]    = ind2sub(size(DesignedVA),find(DesignedVA == 1));

    [yTx, xTx]      = ind2sub(size(tx),find(tx >= 1));
    [yRx, xRx]      = ind2sub(size(rxStar),find(rxStar == 1));

    % Number of Rx
    NrxStar         = sum(nnz(rxStar));

    % Display Results
    fprintf('--------------------------------------------------------------------------------\n');
    fprintf('Number of Tx Elements     = %3d\n', sum(nnz(tx)));
    fprintf('Number of Rx Elements     = %3d\n', NrxStar);
    fprintf('VA Placement Failure Rate = %.1f Percent \n', nnz(DesignedVA - desiredVA.') / nnz(desiredVA) * 100);
    fprintf('================================================================================\n');

    % Plot Results
    figure('Name','Filled Tx-Rx Configuration','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/1.5) floor(ScreenSize(4)/7)]);
    movegui('north')
    plot(xTx,yTx,'o','MarkerSize',5,'LineWidth',1.5), hold on
    plot(xRx,yRx + 1,'v','MarkerSize',5,'LineWidth',1.5)
    xlabel('Horizontal (Half Wavelength)','Interpreter','latex','FontSize',12')
    legend('Tx Chain', 'Rx Chain', 'Interpreter','latex','Location','best','FontSize',11)
    grid on
    xlim([0 max(max(xTx), max(xRx))+1])
    ylim([0.8 2.2])
    yticks([1 2])
    yticklabels({'Tx', 'Rx'})

    figure('Name','Filled Virtual Array Configuration','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/1.5) floor(ScreenSize(4)/7)]);
    movegui('center')
    plot(xVA,yVA,'s','MarkerSize',10,'LineWidth',1.5), hold on
    plot(xdVA,ydVA,'.','MarkerSize',12)
    xlabel('Horizontal (Half Wavelength)','Interpreter','latex','FontSize',12')
    legend('Desired VA', 'Designed VA', 'Interpreter','latex','Location','best','FontSize',11)
    grid on
    xlim([0 max(max(xVA), max(xdVA))+1])
    ylim([0.8 1.2])
    yticks(1)
    yticklabels({'VA'})

    drawnow
    clearvars -except figureNum ScreenSize

end

%%  Figure 17 (Filled VA Design - 100)

if any(17 == figureNum)

    reverseStr      = '';
    lengthVA = 100;
    txPos = E_SMART(50,'result','off');

    % - Antenna Configration
    tx(txPos + 1)   = 1;
    lengthRx        = lengthVA - length(tx) + 1;

    desiredVA       = ones(lengthVA,1);
    minRx           = ceil(lengthVA / length(txPos));

    % - Optimization
    nVar            = lengthRx;
    coeffMat        = double_blocked(tx,[1 lengthRx]);
    epsilon         = 1e-2;
    threshold       = 1e-6;

    subIter         = 0;
    maxIter         = 1e3;
    subIterMax      = 1e2;

    % Pre Allocation
    objfunValP      = zeros(subIterMax,1);
    rxStarMat       = zeros(nVar,maxIter);

    % Main Loop
    fprintf('================================================================================\n');
    fprintf('Optimizing ...\n');
    for iter = 1:maxIter

        subIter = subIter + 1;

        if subIter == 1
            C = eye(nVar);
        else
            C = diag(1./(rx + epsilon));
        end

        cvx_begin quiet

        cvx_solver gurobi

        variable rx(nVar) nonnegative

        minimize sum(pow_abs(coeffMat*rx - desiredVA,2))

        subject to

        ones(1,nVar) * C*rx <= minRx                            %#ok
        rx <= 1 + 1e-3;                                         %#ok
        C * rx <= 1 + 1e-3;                                     %#ok

        cvx_end

        rxStarMat(:,subIter) = C * rx;

        rxOptIter = round(rxStarMat(:,subIter)).';
        designedVAiter = conv2(tx,rxOptIter);
        designedVAiter(designedVAiter >= 1) = 1;

        objfunValP(subIter) = norm(designedVAiter - desiredVA.', 1);

        msg = sprintf('Sub Iteration %2d: Objective Function Value = %3d, Minimum Rx Sensors = %2d \n', ...
            subIter, objfunValP(subIter), minRx);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));

        if norm(designedVAiter - desiredVA.', 1) < 1
            break
        end

        if subIter >= 5 && abs(max(objfunValP(subIter:-1:subIter-4)) - mean(objfunValP(subIter-4:subIter))) <= threshold
            minRx = minRx + 1;
            subIter = 0;
        end

    end
    fprintf('================================================================================\n');
    fprintf('Optimization Completed\n');

    rxStar          = round(rxStarMat(:,subIter)).';

    % - Results
    [yVA, xVA]      = ind2sub(size(desiredVA.'),find(desiredVA.' == 1));

    DesignedVA      = conv2(tx,rxStar);
    DesignedVA(DesignedVA >=1) = 1;
    [ydVA, xdVA]    = ind2sub(size(DesignedVA),find(DesignedVA == 1));

    [yTx, xTx]      = ind2sub(size(tx),find(tx >= 1));
    [yRx, xRx]      = ind2sub(size(rxStar),find(rxStar == 1));

    % Number of Rx
    NrxStar         = sum(nnz(rxStar));

    % Display Results
    fprintf('--------------------------------------------------------------------------------\n');
    fprintf('Number of Tx Elements     = %3d\n', sum(nnz(tx)));
    fprintf('Number of Rx Elements     = %3d\n', NrxStar);
    fprintf('VA Placement Failure Rate = %.1f Percent \n', nnz(DesignedVA - desiredVA.') / nnz(desiredVA) * 100);
    fprintf('================================================================================\n');

    % Plot Results
    figure('Name','Filled Tx-Rx Configuration','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/1.5) floor(ScreenSize(4)/7)]);
    movegui('north')
    plot(xTx,yTx,'o','MarkerSize',5,'LineWidth',1.5), hold on
    plot(xRx,yRx + 1,'v','MarkerSize',5,'LineWidth',1.5)
    xlabel('Horizontal (Half Wavelength)','Interpreter','latex','FontSize',12')
    legend('Tx Chain', 'Rx Chain', 'Interpreter','latex','Location','best','FontSize',11)
    grid on
    xlim([0 max(max(xTx), max(xRx))+1])
    ylim([0.8 2.2])
    yticks([1 2])
    yticklabels({'Tx', 'Rx'})

    figure('Name','Filled Virtual Array Configuration','NumberTitle','off',...
        'Position', [0 0 floor(ScreenSize(3)/1.5) floor(ScreenSize(4)/7)]);
    movegui('center')
    plot(xVA,yVA,'s','MarkerSize',10,'LineWidth',1.5), hold on
    plot(xdVA,ydVA,'.','MarkerSize',12)
    xlabel('Horizontal (Half Wavelength)','Interpreter','latex','FontSize',12')
    legend('Desired VA', 'Designed VA', 'Interpreter','latex','Location','best','FontSize',11)
    grid on
    xlim([0 max(max(xVA), max(xdVA))+1])
    ylim([0.8 1.2])
    yticks(1)
    yticklabels({'VA'})

    drawnow
    clearvars -except figureNum ScreenSize

end

%%  Table 3 & 4 (Wiehgt Function, Coupling, and Leakage)

if any(18 == figureNum)

    load CouplingMat

    N                   = 24;
    DOFsmart            = 156;
    DOFE_SMART          = 149;
    L                   = 16;
    c                   = 3e8;
    fc                  = 77e9;
    lambda              = c / fc;

    % - Array
    % SMART and E-SMART
    indicesSMART        = SMART(DOFsmart,'result','off');
    indicesE_SMART      = E_SMART(DOFE_SMART,'result','off');
    % MISC
    indicesMISC         = MISC(N);
    % ANA-II1
    indicesAug          = AugNestedII1Gen(N);
    % MRA
    indicesMRA          = [0 1 2 3 7 14 21 28 43 58 73 88 103 118 133 148 163 171 179 187 195 196 197 198];
    % CCA
    indicesCCA          = CCA(13,1,12);

    % - Weight Function
    wSMART             = histcounts(pdist(indicesSMART'), 1:max(indicesSMART));
    wE_SMART           = histcounts(pdist(indicesE_SMART'), 1:max(indicesE_SMART));
    wAug               = histcounts(pdist(indicesAug'), 1:max(indicesAug));
    wMISC              = histcounts(pdist(indicesMISC'), 1:max(indicesMISC));
    wMRA               = histcounts(pdist(indicesMRA'), 1:max(indicesMRA));
    wCCA               = histcounts(pdist(indicesCCA'), 1:max(indicesCCA));

    methodName      = {'SMART'
        'E-SMART'
        'ANA-II1'
        'MISC'
        'CCA'
        'MRA'};   

    % Pre Allocation
    Leakage         = zeros(1, length(methodName));    

    for method = 1:length(methodName)

        Leakage(method) = norm(abs(CouplingMat{method}) - diag(diag(abs(CouplingMat{method}))), 'fro') / norm(abs(CouplingMat{method}), 'fro');

        figure('Name', sprintf('%s C Matrix', methodName{method}), ...
            'NumberTitle', 'off', ...
            'Position', [(method-1)*floor(ScreenSize(3)/6) ...
            floor(ScreenSize(4)/2)-floor(ScreenSize(4)/4) ...
            floor(ScreenSize(3)/6) ...
            floor(ScreenSize(3)/6)]);
        imagesc(20*log10(abs(CouplingMat{method})))
        colormap(summer)
        colorbar
        clim([-50 inf])
        axis square
        title(sprintf('%s (Leakage $= %.5f$)', methodName{method}, Leakage(method)), 'Interpreter', 'Latex', 'FontSize', 12);

    end

    % - Results
    figure('Name','w SMART','NumberTitle','off',...
        'Position', [0 1.9*floor(ScreenSize(4)/3) floor(ScreenSize(3)/6) floor(ScreenSize(4)/3.5)]);
    stem(-L:L, [fliplr(wSMART(1:L)) N wSMART(1:L)], 'filled', 'LineWidth',1.5)
    xlabel('Coarray Location $l$','Interpreter','latex','FontSize',15')
    ylabel('$w(l)$','Interpreter','latex','FontSize',15')
    ax = gca;
    ax.FontSize = 14;
    xlim([0 inf])
    xticks(0:4:L)
    ylim([0 N])
    yticks(0:5:N)
    grid on
    box on
    axis square

    figure('Name','w E-SMART','NumberTitle','off',...
        'Position', [floor(ScreenSize(3)/6) 1.9*floor(ScreenSize(4)/3) floor(ScreenSize(3)/6) floor(ScreenSize(4)/3.5)]);
    stem(-L:L, [fliplr(wE_SMART(1:L)) N wE_SMART(1:L)], 'filled', 'LineWidth',1.5)
    xlabel('Coarray Location $l$','Interpreter','latex','FontSize',15')
    ylabel('$w(l)$','Interpreter','latex','FontSize',15')
    ax = gca;
    ax.FontSize = 14;
    xlim([0 inf])
    xticks(0:4:L)
    ylim([0 N])
    yticks(0:5:N)
    grid on
    box on
    axis square

    figure('Name','w ANA-II1','NumberTitle','off',...
        'Position', [2*floor(ScreenSize(3)/6) 1.9*floor(ScreenSize(4)/3) floor(ScreenSize(3)/6) floor(ScreenSize(4)/3.5)]);
    stem(-L:L, [fliplr(wAug(1:L)) N wAug(1:L)], 'filled', 'LineWidth',1.5)
    xlabel('Coarray Location $l$','Interpreter','latex','FontSize',15')
    ylabel('$w(l)$','Interpreter','latex','FontSize',15')
    ax = gca;
    ax.FontSize = 14;
    xlim([0 inf])
    xticks(0:4:L)
    ylim([0 N])
    yticks(0:5:N)
    grid on
    box on
    axis square

    figure('Name','w MISC','NumberTitle','off',...
        'Position', [3*floor(ScreenSize(3)/6) 1.9*floor(ScreenSize(4)/3) floor(ScreenSize(3)/6) floor(ScreenSize(4)/3.5)]);
    stem(-L:L, [fliplr(wMISC(1:L)) N wMISC(1:L)], 'filled', 'LineWidth',1.5)
    xlabel('Coarray Location $l$','Interpreter','latex','FontSize',15')
    ylabel('$w(l)$','Interpreter','latex','FontSize',15')
    ax = gca;
    ax.FontSize = 14;
    xlim([0 inf])
    xticks(0:4:L)
    ylim([0 N])
    yticks(0:5:N)
    grid on
    box on
    axis square

    figure('Name','w CCA','NumberTitle','off',...
        'Position', [4*floor(ScreenSize(3)/6) 1.9*floor(ScreenSize(4)/3) floor(ScreenSize(3)/6) floor(ScreenSize(4)/3.5)]);
    stem(-L:L, [fliplr(wCCA(1:L)) N wCCA(1:L)], 'filled', 'LineWidth',1.5)
    xlabel('Coarray Location $l$','Interpreter','latex','FontSize',15')
    ylabel('$w(l)$','Interpreter','latex','FontSize',15')
    ax = gca;
    ax.FontSize = 14;
    xlim([0 inf])
    xticks(0:4:L)
    ylim([0 N])
    yticks(0:5:N)
    grid on
    box on
    axis square

    figure('Name','w MRA','NumberTitle','off',...
        'Position', [5*floor(ScreenSize(3)/6) 1.9*floor(ScreenSize(4)/3) floor(ScreenSize(3)/6) floor(ScreenSize(4)/3.5)]);
    stem(-L:L, [fliplr(wMRA(1:L)) N wMRA(1:L)], 'filled', 'LineWidth',1.5)
    xlabel('Coarray Location $l$','Interpreter','latex','FontSize',15')
    ylabel('$w(l)$','Interpreter','latex','FontSize',15')
    ax = gca;
    ax.FontSize = 14;
    xlim([0 inf])
    xticks(0:4:L)
    ylim([0 N])
    yticks(0:5:N)
    grid on
    box on
    axis square

    fprintf('===================== Weight =====================\n');
    fprintf('SMART                  : w1 = %2d, w2 = %2d, w3 = %2d\n', wSMART(1), wSMART(2), wSMART(3));
    fprintf('E-SMART                : w1 = %2d, w2 = %2d, w3 = %2d\n', wE_SMART(1), wE_SMART(2), wE_SMART(3));
    fprintf('ANA-II1                : w1 = %2d, w2 = %2d, w3 = %2d\n', wAug(1), wAug(2), wAug(3));
    fprintf('MISC                   : w1 = %2d, w2 = %2d, w3 = %2d\n', wMISC(1), wMISC(2), wMISC(3));
    fprintf('CCA                    : w1 = %2d, w2 = %2d, w3 = %2d\n', wCCA(1), wCCA(2), wCCA(3));
    fprintf('MRA                    : w1 = %2d, w2 = %2d, w3 = %2d\n', wMRA(1), wMRA(2), wMRA(3));
    fprintf('==================================================\n\n');

    fprintf('========== Coupling Leakage ==========\n');
    fprintf('SMART Leakage                = %0.5f\n', Leakage(1));
    fprintf('E-SMART Leakage              = %0.5f\n', Leakage(2));
    fprintf('ANA-II1 Leakage              = %0.5f\n', Leakage(3));
    fprintf('MISC Leakage                 = %0.5f\n', Leakage(4));
    fprintf('CCA Leakage                  = %0.5f\n', Leakage(5));
    fprintf('MRA Leakage                  = %0.5f\n', Leakage(6));
    fprintf('======================================\n');

    drawnow
    clearvars -except figureNum ScreenSize

end