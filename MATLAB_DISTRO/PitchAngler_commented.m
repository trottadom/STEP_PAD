addpath('C:\Users\trotta\Desktop\Shocks'); % Check and then remove
trotta_fuse();                             % Check and then remove
dump_figures = true; % Flaf for figure printing
%% SolO MAG Data
SolO.name = "Solar Orbiter"; % Create SC structure
d_prelim_burst  = 'C:\Users\trotta\OneDrive - Imperial College London\Flight Data\Archive\solo_L2_mag-rtn-normal\2022\'; % Directory where MAG data is stored
addpath(d_prelim_burst);
d_SWA_gm = 'C:\Users\trotta\Desktop\Shocks\Data\SWA\Survey\ground_moments\'; % Directory where SWA data is stored (we don't use it)
d_SWA_ef = "C:\Users\trotta\Desktop\Shocks\Data\SWA\Survey\Eflux\";  % Directory where SWA data is stored (not in use)
d_STEP = '.\data\'; % Directory where STEP dataset is stored
%d_EPT = "C:\Users\trotta\Desktop\Shocks\Data\EPD\Survey\EPT\";
q_day = '*20220415*'; % Day handle
SolO = load_MAG_data(SolO,d_prelim_burst,q_day); % Load MAG data (check SERPENTINE loader if needed https://github.com/serpentine-h2020/serpentine/blob/main/notebooks/sep_analysis_tools/data_loader.ipynb)
% Aim to have time and B arrays
%%
SolO = load_STEP_data(SolO,d_STEP,q_day,"magnet_fluxes"); % Load STEP data (reduntant piece of code, just used for the Energy vector, procedure to be added to load_STEP_elx_pad)
%%
SolO = load_STEP_elx_pad(SolO,d_STEP,q_day); % The actal function we are after (We may want to change the name of this)
%%
% Generalised laziness (can get rid)
SolO.E_pix = SolO.Evec_STEP;
%% Sectors plot (sanity check)
figure(1234);
til = tiledlayout(3,5,'TileSpacing','Compact','Padding','Compact');
xl1 = datetime(2022,04,15,2,30,0); % Start time
xl2 = datetime(2022,04,15,4,30,0); % End time
% Most raw kind of plot
for ipix = 1:15
    nexttile;
    pcolor(squeeze(SolO.t_step),squeeze(SolO.E_pix),(squeeze(SolO.elx_pix_trick(:,:,ipix))')); % x = time_STEP, y = Energy , intensity = flux * E^2
    shading 'flat';
    axx = gca;
    axx.YScale = "log";
    colorbar
    colormap('jet');
    xlim([xl1,xl2])
    clim([2000 52000]);
end
%%
%shock_time_SolO = datetime(2021, 11, 03, 14, 04, 25, 0000, 'Format', 'dd-MMMM-yyyy HH:mm:ss.SSS');
fntsz = 13;
% Select overview stream
yesmom=0;
st_ovi_t  = xl1; % Start time 
end_ovi_t = xl2; % End time
sSolO = select_substream_B_mom(st_ovi_t,end_ovi_t, SolO,yesmom);    % Crop magnetic field data 
sSolO = select_substream_STEP_Elx(st_ovi_t,end_ovi_t, SolO, sSolO); % Crop STEP data
%% Now sync STEP PAS and MAG
eVJ = 1.60218*10^(-19); % electronvolt to Joule
mp  = 1.67262192*10^(-27); % Mass of proton
me  = 9.10938e-31; % Mass of electron

E_STEP = sSolO.E_pix.*10^6; % energy vector in electronvolt  
EJ  = E_STEP.*eVJ;    % Energy in Joule
Vvec = sqrt(EJ.*(2/me)).*10^(-3); % Vector of speeds in km/s
%sSolO.t_STEP = t_step; % generalised laziness
% RESAMPLING (goal: have B vector at same rate as STEP vector)
for it =1:size(sSolO.t_STEP,1)-1
    itB = get_B_indices(sSolO, it);
    itB2  = get_B_indices(sSolO, it+1);
    deltaB = itB2 -itB;
    %deltaM = itm2 -itm;
    sSolO.B_sync(it,:) = mean(sSolO.B(itB:itB+deltaB,:),1);
    %sSolO.U_sync(it,:) = mean(sSolO.U(itm:itm+deltaM,:),1);
end
sSolO.B_sync(size(sSolO.t_STEP,1),:) = sSolO.B_sync(size(sSolO.t_STEP,1)-1,:); % Bsync is the B field at STEP resolution
%sSolO.U_sync(size(sSolO.t_STEP,1),:) = sSolO.U_sync(size(sSolO.t_STEP,1)-1,:);

%% Generate Velocity vectors
% This, at the end, gives you a collection of velocity vectors pointing to
% each pixel and centred on the SC rest frame ()
VRf = zeros(size(sSolO.pixel_RTN,1),size(Vvec,1));
VTf = zeros(size(sSolO.pixel_RTN,1),size(Vvec,1));
VNf = zeros(size(sSolO.pixel_RTN,1),size(Vvec,1));
for iv =1:size(Vvec,1)  % velocity (energy) loop
    for ipix =1:size(sSolO.pixel_RTN,1) % pixel loop
        VRf(ipix,iv) = Vvec(iv)*sSolO.pixel_RTN(ipix,1);
        VTf(ipix,iv) = Vvec(iv)*sSolO.pixel_RTN(ipix,2);
        VNf(ipix,iv) = Vvec(iv)*sSolO.pixel_RTN(ipix,3);
    end
end
%% Now generate PAD using these vectors
for it =1:size(sSolO.t_STEP,1) % time loop
    for iv =1:size(Vvec,1)  % velocity (energy) loop
        for ipix =1:size(sSolO.pixel_RTN,1) % pixel loop
            vR = VRf(ipix,iv); % NOT Corrected for bulk flow speed
            vT = VTf(ipix,iv); 
            vN = VNf(ipix,iv);
            BR = sSolO.B_sync(it,1);  % Magnetic field vector at that time
            BT = sSolO.B_sync(it,2);
            BN = sSolO.B_sync(it,3);
            vm = sqrt(vR^2+vT^2+vN^2);
            Bm = sqrt(BR^2+BT^2+BN^2);
            vdotB = vR*BR + vT*BT + vN*BN; % scalar product
            if abs(vdotB/(vm*Bm)) > 1 % Just for sanity 
                disp("Careful")
            end
            sSolO.pad_f(it,iv,ipix) = acosd(vdotB/(vm*Bm)); % pitch angle distribution Nt x Ech x pixel
        end
    end
end

%% Now try and DO IT FOR ALL ENERGY CHANNELS
figure(112)
tim = 145; % time?
ll = linspace(1,8); % this should span the energy bins, I guess
ech = 2; % Energy channel
szdt1 = 30; % No idea
stp=4; % No idea (one you make sure they dont o much can erase)
Npts = 25; % IMPORTANT Points with which I want to resolve my half-circle
edges = linspace(0,180,Npts);
%edges = 0;
histT = zeros(Npts,1); % Setup histogram
Nhit  = zeros(Npts,1); % I guess how many bins I hit (not sure)
histT2 = zeros(Npts,1); % Setup histogram
Nhit2  = zeros(Npts,1); % I guess how many bins I hit (not sure)

Dalpha = edges(2)-edges(1); % Resolution in pitch angle

h2  = histogram(sSolO.pad_f(tim, ech,:), edges); % This IS CRUCIAL and gives me the coverage of each pixel

% See which pitch angle bins you hit and put in histogram (populate the
% grid)
% Populate the grid
for ipix = 1:15
    ival2 = int32(sSolO.pad_f(tim, ech, ipix)/Dalpha)+1;
    Nhit2(ival2) = Nhit2(ival2)+1;
    histT2(ival2) = histT(ival2) + sSolO.elx_pix_trick(tim,ech,ipix);
end
%
% Simple averaging operation (for cases where more than 1 pixel is in a
% certain grid portion)
for ipt = 1:Npts
    if Nhit2(ipt > 0)
        histT2(ipt) = histT2(ipt)/Nhit2(ipt);
    end
end

% Create mask and timeseries
Nt = size(sSolO.elx_pix_trick,1);
Nch = size(sSolO.elx_pix,2);
sSolO.PAD_GRID = zeros(Nt,Npts,Nch); % This is the grid you are after at the end
sSolO.PAD_GRID(:,:,:) = NaN; 


% Populating the grid for each time AND energy
for ech = 1:Nch
    for it =1:Nt
        histT = zeros(Npts,1);
        Nhit  = zeros(Npts,1);
    
      
        Dalpha = edges(2)-edges(1);
        %h  = histogram(sSolO.pad(tim, ech,:), edges);
        for ipix = 1:15
            ival = int32(sSolO.pad_f(it, ech, ipix)/Dalpha)+1;
            Nhit(ival) = Nhit(ival)+1;
            %histT(ival) = histT(ival) + sSolO.magnet_pix(it,ech,ipix)*sSolO.Evec_STEP(ech)^2;
            histT(ival) = histT(ival) + sSolO.elx_pix(it,ech,ipix);
    
      
        end
        %
        for ipt = 1:Npts
            if Nhit(ipt) > 0
                histT(ipt) = histT(ipt)/Nhit(ipt);
                sSolO.PAD_GRID(it,ipt,ech) = (histT(ipt));
                %disp("ping")
            end
        end
    
     
    
    end
end
disp("Done");



%% Plot

sSolO.PAD_GRID(sSolO.PAD_GRID<=0) = NaN;
SizeX = 780;
SizeY = 900;
FontSize = 13;
screensize   = get(0,'ScreenSize');
ReqDIMENSION = round([SizeX SizeY]);
DIMENSION    = min([ReqDIMENSION;screensize(3:4)]);
FIG.main.axis =figure('Position',[0 (screensize(4)-DIMENSION(2)) DIMENSION]);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize',FontSize)
set(gcf,'color','w');
fntsz = 13;

cmin = 5;  %3.3
cmax = 8;  %3.4
%xl1 = datetime(2020,12,10,7,0,0);
%xl2 = datetime(2020,12,10,9,30,0);


til = tiledlayout(6,1);
ech=1;
nexttile;
    pcolor(sSolO.t_STEP, edges,log10(sSolO.PAD_GRID(:,:,ech))');
    shading 'flat'
    colormap 'jet';
    xlim([xl1, xl2]);
    ylim([0 185]);
    clim([cmin cmax]);
    cbr = colorbar;
    cbr.Label.String = "$\rm log(Diff. \,Flux)$";
    cbr.Label.Interpreter = "latex";
    tt = text(xl1+duration(0,5,0,0), 20, sSolO.pix_ch_txt(ech,:));
    axx = gca;
    axx.XTick = [];
    axx.YLabel.String = "Pitch Angle $[^\circ]$";
    axx.YLabel.Interpreter = "latex";

ech=2;
nexttile;
    pcolor(sSolO.t_STEP, edges,log10(sSolO.PAD_GRID(:,:,ech))');
    shading 'flat'
    colormap 'jet';
    xlim([xl1, xl2]);
    ylim([0 185]);
    clim([cmin cmax]);
    cbr = colorbar;
    cbr.Label.String = "$\rm log(Diff. \, Flux)$";
    cbr.Label.Interpreter = "latex";
    tt = text(xl1+duration(0,5,0,0), 20, sSolO.pix_ch_txt(ech,:));
    axx = gca;
    axx.XTick = [];
    axx.YLabel.String = "Pitch Angle $[^\circ]$";
    axx.YLabel.Interpreter = "latex";
    
ech=3;
nexttile;
    pcolor(sSolO.t_STEP, edges,log10(sSolO.PAD_GRID(:,:,ech))');
    shading 'flat'
    colormap 'jet';
    xlim([xl1, xl2]);
    ylim([0 185]);
    clim([cmin cmax]);
    cbr = colorbar;
    cbr.Label.String = "$\rm log(Diff. \, Flux)$";
    cbr.Label.Interpreter = "latex";
    tt = text(xl1+duration(0,5,0,0), 20, sSolO.pix_ch_txt(ech,:));
    axx = gca;
    axx.XTick = [];
    axx.YLabel.String = "Pitch Angle $[^\circ]$";
    axx.YLabel.Interpreter = "latex";

ech=4;
nexttile;
    pcolor(sSolO.t_STEP, edges,log10(sSolO.PAD_GRID(:,:,ech))');
    shading 'flat'
    colormap 'jet';
    xlim([xl1, xl2]);
    ylim([0 185]);
    clim([cmin cmax]);
    cbr = colorbar;
    cbr.Label.String = "$\rm log(Diff. \, Flux)$";
    cbr.Label.Interpreter = "latex";
    tt = text(xl1+duration(0,5,0,0), 20, sSolO.pix_ch_txt(ech,:));
    axx = gca;
    axx.XTick = [];
    axx.YLabel.String = "Pitch Angle $[^\circ]$";
    axx.YLabel.Interpreter = "latex";

ech=10;
nexttile;
    pcolor(sSolO.t_STEP, edges,log10(sSolO.PAD_GRID(:,:,ech))');
    shading 'flat'
    colormap 'jet';
    xlim([xl1, xl2]);
    ylim([0 185]);
    clim([cmin cmax]);
    cbr = colorbar;
    cbr.Label.String = "$\rm log(Diff. \, Flux)$";
    cbr.Label.Interpreter = "latex";
    tt = text(xl1+duration(0,5,0,0), 20, sSolO.pix_ch_txt(ech,:));
    axx = gca;
    axx.XTick = [];
    axx.YLabel.String = "Pitch Angle $[^\circ]$";
    axx.YLabel.Interpreter = "latex";

%ech=5;
nexttile;
    p1 = plot(sSolO.tB,sSolO.B(:,1),'Color', [.3, .3, .8],'LineWidth',1.2);
    hold on;
    p2 = plot(sSolO.tB,sSolO.B(:,2),'Color', [.3, .8, .3],'LineWidth',1.2);
    p3 = plot(sSolO.tB,sSolO.B(:,3),'Color', [.8, .3, .3],'LineWidth',1.2);
    p4 = plot(sSolO.tB,sSolO.B(:,4),'Color', [0 0 0],'LineWidth',1.2);
    hold off;
    xlim([xl1, xl2]);
    ylim([-12 13]);
    yl = yline(0);
    yl.Color = [.1 .1 .1];
    yl.LineWidth = 1.3;
    yl.LineStyle = '-.';
    %clim([cmin cmax]);
    %cbr = colorbar;
    %cbr.Label.String = "$\rm log(Diff. \, Flux)$";
    %cbr.Label.Interpreter = "latex";
    %tt = text(xl1+duration(0,5,0,0), 20, pix_ch_txt(ech,:));
    ll = legend([p1 p2 p3 p4], "$\rm B_R$", "$\rm B_T$", "$\rm B_N$", "$\rm B$");
    ll.Interpreter = "latex";
    ll.Location = "northwest";
    axx = gca;
    axx.XLabel.String = "UT";
    axx.XLabel.Interpreter = "latex";

    axx.YLabel.String = "B [nT]";
    axx.YLabel.Interpreter = "latex";


%print(gcf,"Camille_10Dec_Prelim_PADs.jpg",'-djpeg');


















